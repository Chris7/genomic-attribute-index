use crate::*;

#[path = "query/bed.rs"]
mod bed;
#[path = "query/gff.rs"]
mod gff;

enum IndexStorage {
    Owned(Vec<u8>),
    Mapped(Mmap),
}

impl AsRef<[u8]> for IndexStorage {
    fn as_ref(&self) -> &[u8] {
        match self {
            Self::Owned(bytes) => bytes,
            Self::Mapped(bytes) => bytes,
        }
    }
}

/// A bounds-checked GAI reader. It only decompresses the posting and span
/// blocks needed by an individual lookup. [`Self::open_mmap`] keeps the FST
/// and fixed directories in an OS memory map instead of copying the file.
pub struct NameIndexReader {
    bytes: IndexStorage,
    metadata: IndexMetadata,
    sections: BTreeMap<SectionKind, SectionDirectoryEntry>,
    posting_directory: Vec<PostingDirectoryEntry>,
    span_directory: Vec<SpanDirectoryEntry>,
}

impl NameIndexReader {
    /// Opens and validates a GAI file without opening its source files.
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        let bytes = fs::read(path)?;
        Self::from_bytes(bytes)
    }

    /// Opens and validates a GAI using a read-only memory map.
    pub fn open_mmap(path: impl AsRef<Path>) -> Result<Self> {
        let file = File::open(path)?;
        // SAFETY: the file descriptor remains valid while the map is created;
        // the resulting Mmap owns the mapping for the reader's lifetime.
        let bytes = unsafe { Mmap::map(&file)? };
        Self::from_storage(IndexStorage::Mapped(bytes))
    }

    /// Parses a GAI byte stream.  This is useful for corruption tests and for
    /// callers that memory-map the file themselves before handing it to GAI.
    pub fn from_bytes(bytes: Vec<u8>) -> Result<Self> {
        Self::from_storage(IndexStorage::Owned(bytes))
    }

    fn from_storage(storage: IndexStorage) -> Result<Self> {
        let bytes = storage.as_ref();
        if bytes.len() < HEADER_SIZE {
            return Err(Error::Corrupt("truncated GAI header".into()));
        }
        if bytes[..4] != MAGIC {
            return Err(Error::Corrupt("invalid GAI magic".into()));
        }
        let mut offset = 4;
        let major_version = read_u16(bytes, &mut offset, "major version")?;
        let minor_version = read_u16(bytes, &mut offset, "minor version")?;
        if major_version != MAJOR_VERSION {
            return Err(Error::Corrupt(format!(
                "unsupported GAI major version {major_version}"
            )));
        }
        if minor_version > MINOR_VERSION {
            return Err(Error::Corrupt(format!(
                "unsupported GAI minor version {minor_version}"
            )));
        }
        let flags = read_u32(bytes, &mut offset, "flags")?;
        if flags & !1 != 0 {
            return Err(Error::Corrupt("unknown GAI flags".into()));
        }
        let byte_order = read_u8(bytes, &mut offset, "byte order")?;
        let coordinate_convention = read_u8(bytes, &mut offset, "coordinate convention")?;
        let normalization = read_u8(bytes, &mut offset, "normalization policy")?;
        let reserved = read_u8(bytes, &mut offset, "header reserved byte")?;
        if byte_order != BYTE_ORDER_LITTLE
            || coordinate_convention != COORDINATE_ZERO_BASED_HALF_OPEN
            || reserved != 0
        {
            return Err(Error::Corrupt("unsupported GAI header conventions".into()));
        }
        if normalization > NORMALIZATION_CASE_SENSITIVE {
            return Err(Error::Corrupt("unknown normalization policy".into()));
        }
        let header_size = read_u32(bytes, &mut offset, "header size")? as usize;
        let directory_entry_size = read_u32(bytes, &mut offset, "directory entry size")? as usize;
        let section_count = read_u32(bytes, &mut offset, "section count")? as usize;
        let reserved = read_u32(bytes, &mut offset, "header reserved field")?;
        if header_size != HEADER_SIZE
            || directory_entry_size != DIRECTORY_ENTRY_SIZE
            || reserved != 0
        {
            return Err(Error::Corrupt("unsupported GAI header layout".into()));
        }
        if section_count != 7 {
            return Err(Error::Corrupt("unsupported GAI section count".into()));
        }
        let term_count = read_u64(bytes, &mut offset, "term count")?;
        let unique_span_count = read_u64(bytes, &mut offset, "span count")?;
        let posting_count = read_u64(bytes, &mut offset, "posting count")?;
        let postings_block_count = read_u64(bytes, &mut offset, "postings block count")?;
        let span_block_count = read_u64(bytes, &mut offset, "span block count")?;
        let max_posting_count = term_count
            .checked_mul(unique_span_count)
            .ok_or_else(|| Error::Corrupt("posting count overflow".into()))?;
        if posting_count > max_posting_count
            || (term_count == 0 && posting_count != 0)
            || (term_count > 0 && posting_count < term_count)
            || (term_count > 0 && postings_block_count == 0)
            || (unique_span_count > 0 && span_block_count == 0)
        {
            return Err(Error::Corrupt("inconsistent GAI counts".into()));
        }
        let source_fingerprint = read_array::<32>(bytes, &mut offset, "source fingerprint")?;
        let coordinate_index_fingerprint =
            read_array::<32>(bytes, &mut offset, "coordinate-index fingerprint")?;
        let reference_dictionary_fingerprint =
            read_array::<32>(bytes, &mut offset, "reference-dictionary fingerprint")?;
        let attribute_count = read_u32(bytes, &mut offset, "attribute count")? as usize;
        let spans_per_block = read_u32(bytes, &mut offset, "span block size")?;
        if spans_per_block == 0 || spans_per_block > 1_000_000 {
            return Err(Error::Corrupt("invalid span block size".into()));
        }
        let section_directory_offset = read_u64(bytes, &mut offset, "section directory offset")?;
        let section_directory_length = read_u64(bytes, &mut offset, "section directory length")?;
        let file_size = read_u64(bytes, &mut offset, "file size")?;
        let reference_count = read_u32(bytes, &mut offset, "reference count")?;
        if reference_count > 1_000_000 {
            return Err(Error::Corrupt("excessive reference count".into()));
        }
        if bytes[204..HEADER_SIZE].iter().any(|byte| *byte != 0) {
            return Err(Error::Corrupt("nonzero reserved header bytes".into()));
        }
        if file_size != bytes.len() as u64 {
            return Err(Error::Corrupt("GAI file size mismatch".into()));
        }
        let directory_end = section_directory_offset
            .checked_add(section_directory_length)
            .ok_or_else(|| Error::Corrupt("section directory overflow".into()))?;
        if section_directory_offset < HEADER_SIZE as u64
            || directory_end > bytes.len() as u64
            || section_directory_length
                != u64::from(section_count as u32)
                    .checked_mul(DIRECTORY_ENTRY_SIZE as u64)
                    .ok_or_else(|| Error::Corrupt("section directory size overflow".into()))?
        {
            return Err(Error::Corrupt("invalid section directory range".into()));
        }
        let mut sections = BTreeMap::new();
        let mut section_ranges = Vec::with_capacity(section_count);
        let mut directory_offset = section_directory_offset as usize;
        for _ in 0..section_count {
            let kind =
                SectionKind::try_from(read_u32(bytes, &mut directory_offset, "section kind")?)?;
            let flags = read_u32(bytes, &mut directory_offset, "section flags")?;
            if flags != 0 {
                return Err(Error::Corrupt("unknown section flags".into()));
            }
            let section_offset = read_u64(bytes, &mut directory_offset, "section offset")?;
            let section_length = read_u64(bytes, &mut directory_offset, "section length")?;
            let item_count = read_u64(bytes, &mut directory_offset, "section item count")?;
            let section_checksum = read_u32(bytes, &mut directory_offset, "section checksum")?;
            let reserved = read_u32(bytes, &mut directory_offset, "section reserved")?;
            if reserved != 0 || section_length > MAX_SECTION_BYTES {
                return Err(Error::Corrupt("invalid section metadata".into()));
            }
            let section_end = section_offset
                .checked_add(section_length)
                .ok_or_else(|| Error::Corrupt("section range overflow".into()))?;
            if section_offset < directory_end || section_end > bytes.len() as u64 {
                return Err(Error::Corrupt("section is outside file".into()));
            }
            section_ranges.push((section_offset, section_end));
            let section = bytes
                .get(section_offset as usize..section_end as usize)
                .ok_or_else(|| Error::Corrupt("section range is not addressable".into()))?;
            if checksum(section) != section_checksum {
                return Err(Error::Corrupt(format!(
                    "{kind:?} section checksum mismatch"
                )));
            }
            if sections
                .insert(
                    kind,
                    SectionDirectoryEntry {
                        kind,
                        flags,
                        offset: section_offset,
                        length: section_length,
                        item_count,
                        checksum: section_checksum,
                    },
                )
                .is_some()
            {
                return Err(Error::Corrupt("duplicate section directory entry".into()));
            }
        }
        section_ranges.sort_unstable();
        if section_ranges
            .windows(2)
            .any(|ranges| ranges[0].1 > ranges[1].0)
        {
            return Err(Error::Corrupt("overlapping GAI sections".into()));
        }
        let required = [
            SectionKind::Attributes,
            SectionKind::Terms,
            SectionKind::PostingsDirectory,
            SectionKind::PostingsData,
            SectionKind::SpansDirectory,
            SectionKind::StartsData,
            SectionKind::LengthsData,
        ];
        if required.iter().any(|kind| !sections.contains_key(kind)) {
            return Err(Error::Corrupt("missing required GAI section".into()));
        }
        let postings_data_length = sections
            .get(&SectionKind::PostingsData)
            .map(|section| section.length)
            .ok_or_else(|| Error::Corrupt("missing postings data section".into()))?;
        let starts_data_length = sections
            .get(&SectionKind::StartsData)
            .map(|section| section.length)
            .ok_or_else(|| Error::Corrupt("missing starts data section".into()))?;
        let lengths_data_length = sections
            .get(&SectionKind::LengthsData)
            .map(|section| section.length)
            .ok_or_else(|| Error::Corrupt("missing lengths data section".into()))?;
        let expected_items = [
            (SectionKind::Attributes, attribute_count as u64),
            (SectionKind::Terms, term_count),
            (SectionKind::PostingsDirectory, postings_block_count),
            (SectionKind::PostingsData, postings_data_length),
            (SectionKind::SpansDirectory, span_block_count),
            (SectionKind::StartsData, starts_data_length),
            (SectionKind::LengthsData, lengths_data_length),
        ];
        for (kind, expected) in expected_items {
            if sections
                .get(&kind)
                .is_some_and(|section| section.item_count != expected)
            {
                return Err(Error::Corrupt(format!("{kind:?} item count mismatch")));
            }
        }
        let section_bytes = |kind: SectionKind| -> Result<&[u8]> {
            let section = sections
                .get(&kind)
                .ok_or_else(|| Error::Corrupt("missing required section".into()))?;
            let end = section
                .offset
                .checked_add(section.length)
                .ok_or_else(|| Error::Corrupt("section range overflow".into()))?;
            bytes
                .get(section.offset as usize..end as usize)
                .ok_or_else(|| Error::Corrupt("section range out of bounds".into()))
        };
        let attributes = decode_attributes(section_bytes(SectionKind::Attributes)?)?;
        if attributes.len() != attribute_count || attributes.is_empty() {
            return Err(Error::Corrupt("attribute count mismatch".into()));
        }
        let terms_bytes = section_bytes(SectionKind::Terms)?;
        let term_map = fst::Map::new(terms_bytes)
            .map_err(|error| Error::Corrupt(format!("invalid term FST: {error}")))?;
        if term_map.len() as u64 != term_count {
            return Err(Error::Corrupt("term count mismatch".into()));
        }
        let posting_directory = decode_posting_directory(
            section_bytes(SectionKind::PostingsDirectory)?,
            section_bytes(SectionKind::PostingsData)?.len() as u64,
            postings_block_count,
        )?;
        let span_directory = decode_span_directory(
            section_bytes(SectionKind::SpansDirectory)?,
            section_bytes(SectionKind::StartsData)?.len() as u64,
            section_bytes(SectionKind::LengthsData)?.len() as u64,
            span_block_count,
            unique_span_count,
            spans_per_block,
            reference_count,
        )?;
        let section_length =
            |kind: SectionKind| -> u64 { sections.get(&kind).map_or(0, |section| section.length) };
        let postings_uncompressed_bytes = posting_directory
            .iter()
            .try_fold(0_u64, |total, entry| {
                total.checked_add(u64::from(entry.uncompressed_length))
            })
            .ok_or_else(|| Error::Corrupt("postings uncompressed size overflow".into()))?;
        let starts_uncompressed_bytes = span_directory
            .iter()
            .try_fold(0_u64, |total, entry| {
                total.checked_add(u64::from(entry.starts_uncompressed_length))
            })
            .ok_or_else(|| Error::Corrupt("starts uncompressed size overflow".into()))?;
        let lengths_uncompressed_bytes = span_directory
            .iter()
            .try_fold(0_u64, |total, entry| {
                total.checked_add(u64::from(entry.lengths_uncompressed_length))
            })
            .ok_or_else(|| Error::Corrupt("lengths uncompressed size overflow".into()))?;
        let span_uncompressed_bytes = starts_uncompressed_bytes
            .checked_add(lengths_uncompressed_bytes)
            .ok_or_else(|| Error::Corrupt("span uncompressed size overflow".into()))?;
        let starts_compressed_blocks = span_directory
            .iter()
            .filter(|entry| entry.starts_compression == 1)
            .count() as u64;
        let lengths_compressed_blocks = span_directory
            .iter()
            .filter(|entry| entry.lengths_compression == 1)
            .count() as u64;
        let file_size = bytes.len() as u64;
        Ok(Self {
            bytes: storage,
            metadata: IndexMetadata {
                major_version,
                minor_version,
                case_sensitive: normalization == NORMALIZATION_CASE_SENSITIVE,
                attributes,
                source_fingerprint,
                coordinate_index_fingerprint,
                reference_dictionary_fingerprint,
                term_count,
                unique_span_count,
                posting_count,
                postings_block_count,
                span_block_count,
                reference_count,
                span_block_size: spans_per_block,
                file_size,
                attribute_section_bytes: section_length(SectionKind::Attributes),
                term_dictionary_bytes: section_length(SectionKind::Terms),
                postings_directory_bytes: section_length(SectionKind::PostingsDirectory),
                postings_uncompressed_bytes,
                postings_data_bytes: section_length(SectionKind::PostingsData),
                span_directory_bytes: section_length(SectionKind::SpansDirectory),
                span_uncompressed_bytes,
                starts_data_bytes: section_length(SectionKind::StartsData),
                lengths_data_bytes: section_length(SectionKind::LengthsData),
                starts_uncompressed_bytes,
                lengths_uncompressed_bytes,
                compressed_postings_blocks: posting_directory
                    .iter()
                    .filter(|entry| entry.compression == 1)
                    .count() as u64,
                delta_start_blocks: span_directory.len() as u64,
                varint_length_blocks: span_directory
                    .iter()
                    .filter(|entry| entry.length_encoding == 0)
                    .count() as u64,
                for_length_blocks: span_directory
                    .iter()
                    .filter(|entry| entry.length_encoding == 2)
                    .count() as u64,
                compressed_start_blocks: starts_compressed_blocks,
                compressed_length_blocks: lengths_compressed_blocks,
            },
            sections,
            posting_directory,
            span_directory,
        })
    }

    /// Returns format metadata and configured attributes.
    pub fn metadata(&self) -> &IndexMetadata {
        &self.metadata
    }

    /// Returns a copy of metadata for inspection interfaces.
    pub fn inspect(&self) -> IndexMetadata {
        self.metadata.clone()
    }

    /// Looks up a normalized term and returns its exact spans.
    pub fn lookup_spans(&self, term: &str) -> Result<Vec<Span>> {
        self.lookup_spans_with_mode(term, MatchMode::Exact)
    }

    /// Looks up a normalized term using the requested match mode and returns
    /// its coordinate-ordered spans. Prefix matching uses the FST automaton;
    /// contains and regex matching scan keys incrementally. Each referenced
    /// postings block is decoded at most once.
    pub fn lookup_spans_with_mode(&self, term: &str, mode: MatchMode) -> Result<Vec<Span>> {
        let (span_ids, _) = self.lookup_span_ids_with_mode_and_stats(term, mode)?;
        self.resolve_span_ids(&span_ids)
    }

    /// Looks up spans and reports the independently decoded bytes used by an
    /// exact query.
    pub fn lookup_spans_with_stats(&self, term: &str) -> Result<(Vec<Span>, LookupStats)> {
        self.lookup_spans_with_mode_and_stats(term, MatchMode::Exact)
    }

    /// Looks up spans with an explicit match mode and reports the independently
    /// decoded bytes used by the lookup. Each postings and span block is
    /// counted once even when several matching terms or span IDs share it.
    pub fn lookup_spans_with_mode_and_stats(
        &self,
        term: &str,
        mode: MatchMode,
    ) -> Result<(Vec<Span>, LookupStats)> {
        let (span_ids, postings_bytes_decompressed) =
            self.lookup_span_ids_with_mode_and_stats(term, mode)?;
        let (spans, _, span_bytes_decompressed) = self.resolve_span_ids_with_stats(&span_ids)?;
        Ok((
            spans,
            LookupStats {
                postings_bytes_decompressed,
                span_bytes_decompressed,
            },
        ))
    }

    /// Returns the sorted span IDs for an exact normalized term.
    pub fn lookup_span_ids(&self, term: &str) -> Result<Vec<u64>> {
        self.lookup_span_ids_with_mode(term, MatchMode::Exact)
    }

    /// Returns the sorted, deduplicated span IDs for a normalized term and
    /// explicit match mode.
    pub fn lookup_span_ids_with_mode(&self, term: &str, mode: MatchMode) -> Result<Vec<u64>> {
        self.lookup_span_ids_with_mode_and_stats(term, mode)
            .map(|(span_ids, _)| span_ids)
    }

    fn lookup_span_ids_with_mode_and_stats(
        &self,
        term: &str,
        mode: MatchMode,
    ) -> Result<(Vec<u64>, u64)> {
        let matcher = CompiledMatch::new(term, mode, self.metadata.case_sensitive)?;
        self.lookup_span_ids_with_matcher_and_stats(&matcher)
    }

    fn lookup_span_ids_with_matcher_and_stats(
        &self,
        matcher: &CompiledMatch,
    ) -> Result<(Vec<u64>, u64)> {
        if matcher.is_empty() {
            return Ok((Vec::new(), 0));
        }
        let terms_section = self.section(SectionKind::Terms)?;
        let term_map = fst::Map::new(terms_section)
            .map_err(|error| Error::Corrupt(format!("invalid term FST: {error}")))?;
        let locators = match matcher {
            CompiledMatch::Exact(query) => term_map.get(query).into_iter().collect(),
            CompiledMatch::Prefix(query) => fst::IntoStreamer::into_stream(
                term_map.search(fst::automaton::Str::new(query).starts_with()),
            )
            .into_values(),
            CompiledMatch::Contains(_) | CompiledMatch::Regex { .. } => {
                let mut locators = Vec::new();
                let mut stream = term_map.stream();
                while let Some((key, locator)) = stream.next() {
                    let key = std::str::from_utf8(key)
                        .map_err(|_| Error::Corrupt("term FST key is not UTF-8".into()))?;
                    if matcher.matches_normalized(key) {
                        locators.push(locator);
                    }
                }
                locators
            }
        };
        if locators.is_empty() {
            return Ok((Vec::new(), 0));
        }

        let mut posting_blocks = BTreeMap::<usize, Vec<u8>>::new();
        let mut span_ids = HashSet::new();
        let mut postings_bytes_decompressed = 0_u64;
        for locator in locators {
            let block_id = usize::try_from(locator >> 32)
                .map_err(|_| Error::Corrupt("posting block ID overflows usize".into()))?;
            let record_offset = usize::try_from(locator & u64::from(u32::MAX))
                .map_err(|_| Error::Corrupt("posting offset overflows usize".into()))?;
            if let std::collections::btree_map::Entry::Vacant(entry) =
                posting_blocks.entry(block_id)
            {
                let bytes = self.decode_posting_block(block_id)?;
                postings_bytes_decompressed = postings_bytes_decompressed
                    .checked_add(bytes.len() as u64)
                    .ok_or(Error::InvalidCoordinate)?;
                entry.insert(bytes);
            }
            let posting_bytes = posting_blocks
                .get(&block_id)
                .ok_or_else(|| Error::Corrupt("missing decoded postings block".into()))?;
            span_ids.extend(decode_posting_record(
                posting_bytes,
                record_offset,
                self.metadata.unique_span_count,
            )?);
        }
        let mut span_ids = span_ids.into_iter().collect::<Vec<_>>();
        span_ids.sort_unstable();
        Ok((span_ids, postings_bytes_decompressed))
    }

    /// Resolves one span ID using a binary search over fixed-width block entries.
    pub fn resolve_span_id(&self, span_id: u64) -> Result<Span> {
        self.resolve_span_ids(&[span_id])?
            .into_iter()
            .next()
            .ok_or_else(|| Error::Corrupt("missing span ID".into()))
    }

    /// Resolves a batch of span IDs while decoding each referenced span block
    /// at most once. The input and output retain the same order.
    pub fn resolve_span_ids(&self, span_ids: &[u64]) -> Result<Vec<Span>> {
        self.resolve_span_ids_with_stats(span_ids)
            .map(|(spans, _, _)| spans)
    }

    fn resolve_span_ids_with_stats(&self, span_ids: &[u64]) -> Result<(Vec<Span>, u64, u64)> {
        let mut requests = BTreeMap::<usize, Vec<(usize, usize)>>::new();
        for (output_index, &span_id) in span_ids.iter().enumerate() {
            if span_id >= self.metadata.unique_span_count {
                return Err(Error::Corrupt("invalid span ID".into()));
            }
            let block_index = self
                .span_directory
                .partition_point(|entry| entry.first_span_id <= span_id)
                .checked_sub(1)
                .ok_or_else(|| Error::Corrupt("span ID is outside span directory".into()))?;
            let entry = &self.span_directory[block_index];
            let row = usize::try_from(span_id - entry.first_span_id)
                .map_err(|_| Error::Corrupt("span row overflows usize".into()))?;
            if row >= entry.span_count as usize {
                return Err(Error::Corrupt("span ID is outside span block".into()));
            }
            requests
                .entry(block_index)
                .or_default()
                .push((output_index, row));
        }

        let mut spans = vec![None; span_ids.len()];
        let mut span_bytes_decompressed = 0_u64;
        for (block_index, rows) in requests.iter() {
            let entry = &self.span_directory[*block_index];
            let starts_payload = self.decode_starts_block(*block_index)?;
            let lengths_payload = self.decode_lengths_block(*block_index)?;
            span_bytes_decompressed = span_bytes_decompressed
                .checked_add(starts_payload.len() as u64)
                .and_then(|total| total.checked_add(lengths_payload.len() as u64))
                .ok_or(Error::InvalidCoordinate)?;
            let decoded = decode_span_rows(&starts_payload, &lengths_payload, entry)?;
            for &(output_index, row) in rows {
                let span = *decoded
                    .get(row)
                    .ok_or_else(|| Error::Corrupt("span row is out of bounds".into()))?;
                spans[output_index] = Some(span);
            }
        }

        let spans = spans
            .into_iter()
            .map(|span| span.ok_or_else(|| Error::Corrupt("missing span ID".into())))
            .collect::<Result<Vec<_>>>()?;
        Ok((spans, requests.len() as u64, span_bytes_decompressed))
    }

    fn section(&self, kind: SectionKind) -> Result<&[u8]> {
        let section = self
            .sections
            .get(&kind)
            .ok_or_else(|| Error::Corrupt("missing GAI section".into()))?;
        let end = section
            .offset
            .checked_add(section.length)
            .ok_or_else(|| Error::Corrupt("section range overflow".into()))?;
        self.bytes
            .as_ref()
            .get(section.offset as usize..end as usize)
            .ok_or_else(|| Error::Corrupt("section range out of bounds".into()))
    }

    fn decode_posting_block(&self, block_id: usize) -> Result<Vec<u8>> {
        let entry = self
            .posting_directory
            .get(block_id)
            .ok_or_else(|| Error::Corrupt("invalid posting block ID".into()))?;
        let data = self.section(SectionKind::PostingsData)?;
        let end = entry
            .compressed_offset
            .checked_add(u64::from(entry.compressed_length))
            .ok_or_else(|| Error::Corrupt("posting block range overflow".into()))?;
        let compressed = data
            .get(entry.compressed_offset as usize..end as usize)
            .ok_or_else(|| Error::Corrupt("posting block range out of bounds".into()))?;
        decompress_block(
            compressed,
            entry.compression,
            entry.uncompressed_length,
            entry.checksum,
            "posting block",
        )
    }

    fn decode_starts_block(&self, block_id: usize) -> Result<Vec<u8>> {
        let entry = self
            .span_directory
            .get(block_id)
            .ok_or_else(|| Error::Corrupt("invalid span block ID".into()))?;
        let data = self.section(SectionKind::StartsData)?;
        let end = entry
            .starts_compressed_offset
            .checked_add(u64::from(entry.starts_compressed_length))
            .ok_or_else(|| Error::Corrupt("span block range overflow".into()))?;
        let compressed = data
            .get(entry.starts_compressed_offset as usize..end as usize)
            .ok_or_else(|| Error::Corrupt("starts block range out of bounds".into()))?;
        decompress_block(
            compressed,
            entry.starts_compression,
            entry.starts_uncompressed_length,
            entry.starts_checksum,
            "span starts block",
        )
    }

    fn decode_lengths_block(&self, block_id: usize) -> Result<Vec<u8>> {
        let entry = self
            .span_directory
            .get(block_id)
            .ok_or_else(|| Error::Corrupt("invalid span block ID".into()))?;
        let data = self.section(SectionKind::LengthsData)?;
        let end = entry
            .lengths_compressed_offset
            .checked_add(u64::from(entry.lengths_compressed_length))
            .ok_or_else(|| Error::Corrupt("span lengths range overflow".into()))?;
        let compressed = data
            .get(entry.lengths_compressed_offset as usize..end as usize)
            .ok_or_else(|| Error::Corrupt("lengths block range out of bounds".into()))?;
        decompress_block(
            compressed,
            entry.lengths_compression,
            entry.lengths_uncompressed_length,
            entry.lengths_checksum,
            "span lengths block",
        )
    }
}

pub(crate) fn decode_posting_record(
    bytes: &[u8],
    record_offset: usize,
    span_count_limit: u64,
) -> Result<Vec<u64>> {
    let mut offset = record_offset;
    let count = read_varint(bytes, &mut offset, "posting span count")?;
    if count > span_count_limit || count > 100_000_000 || count > MAX_BLOCK_BYTES / 8 {
        return Err(Error::Corrupt("excessive posting span count".into()));
    }
    decode_delta_posting_record(bytes, &mut offset, count, span_count_limit)
}

fn decode_delta_posting_record(
    bytes: &[u8],
    offset: &mut usize,
    count: u64,
    span_count_limit: u64,
) -> Result<Vec<u64>> {
    let count = usize::try_from(count)
        .map_err(|_| Error::Corrupt("posting span count overflows usize".into()))?;
    let mut ids = Vec::with_capacity(count);
    let mut previous = 0_u64;
    for index in 0..count {
        let value = read_varint(bytes, offset, "posting span ID")?;
        let id = if index == 0 {
            value
        } else {
            if value == 0 {
                return Err(Error::Corrupt("posting span IDs are not sorted".into()));
            }
            previous
                .checked_add(value)
                .ok_or_else(|| Error::Corrupt("posting span ID overflow".into()))?
        };
        if id >= span_count_limit {
            return Err(Error::Corrupt("posting references invalid span ID".into()));
        }
        ids.push(id);
        previous = id;
    }
    Ok(ids)
}

fn decode_posting_directory(
    bytes: &[u8],
    data_length: u64,
    expected_count: u64,
) -> Result<Vec<PostingDirectoryEntry>> {
    if expected_count > 100_000_000 {
        return Err(Error::Corrupt("excessive postings block count".into()));
    }
    if !bytes.len().is_multiple_of(POSTINGS_DIRECTORY_ENTRY_SIZE)
        || bytes.len() / POSTINGS_DIRECTORY_ENTRY_SIZE != expected_count as usize
    {
        return Err(Error::Corrupt("postings directory size mismatch".into()));
    }
    let mut offset = 0;
    let mut entries = Vec::with_capacity(expected_count as usize);
    let mut previous_end = 0_u64;
    for _ in 0..expected_count {
        let compressed_offset = read_u64(bytes, &mut offset, "posting compressed offset")?;
        let compressed_length = read_u32(bytes, &mut offset, "posting compressed length")?;
        let uncompressed_length = read_u32(bytes, &mut offset, "posting uncompressed length")?;
        let checksum = read_u32(bytes, &mut offset, "posting checksum")?;
        let compression = read_u8(bytes, &mut offset, "posting compression")?;
        let reserved = read_array::<3>(bytes, &mut offset, "posting reserved")?;
        let reserved_tail = read_u64(bytes, &mut offset, "posting reserved tail")?;
        if reserved != [0; 3] || reserved_tail != 0 || compression > 1 {
            return Err(Error::Corrupt("invalid postings directory entry".into()));
        }
        if u64::from(uncompressed_length) > MAX_BLOCK_BYTES {
            return Err(Error::Corrupt("postings block is too large".into()));
        }
        let end = compressed_offset
            .checked_add(u64::from(compressed_length))
            .ok_or_else(|| Error::Corrupt("postings block range overflow".into()))?;
        if compressed_offset < previous_end || end > data_length {
            return Err(Error::Corrupt(
                "postings block is outside data section".into(),
            ));
        }
        previous_end = end;
        entries.push(PostingDirectoryEntry {
            compressed_offset,
            compressed_length,
            uncompressed_length,
            checksum,
            compression,
        });
    }
    Ok(entries)
}

fn decode_span_directory(
    bytes: &[u8],
    starts_data_length: u64,
    lengths_data_length: u64,
    expected_count: u64,
    expected_span_count: u64,
    spans_per_block: u32,
    reference_count: u32,
) -> Result<Vec<SpanDirectoryEntry>> {
    if expected_count > 100_000_000 || expected_span_count > 100_000_000 {
        return Err(Error::Corrupt("excessive span directory count".into()));
    }
    if !bytes.len().is_multiple_of(SPAN_DIRECTORY_ENTRY_SIZE)
        || bytes.len() / SPAN_DIRECTORY_ENTRY_SIZE != expected_count as usize
    {
        return Err(Error::Corrupt("span directory size mismatch".into()));
    }
    let mut offset = 0;
    let mut entries = Vec::with_capacity(expected_count as usize);
    let mut previous_span_id = 0_u64;
    let mut previous_starts_data_end = 0_u64;
    let mut previous_lengths_data_end = 0_u64;
    let mut previous_reference = None;
    let mut previous_start = 0_u64;
    let mut total_rows = 0_u64;
    for index in 0..expected_count {
        let first_span_id = read_u64(bytes, &mut offset, "span first ID")?;
        let span_count = read_u32(bytes, &mut offset, "span count")?;
        let reference_id = read_u32(bytes, &mut offset, "span reference ID")?;
        let first_start = read_u64(bytes, &mut offset, "span first start")?;
        let starts_compressed_offset =
            read_u64(bytes, &mut offset, "span starts compressed offset")?;
        let starts_compressed_length =
            read_u32(bytes, &mut offset, "span starts compressed length")?;
        let starts_uncompressed_length =
            read_u32(bytes, &mut offset, "span starts uncompressed length")?;
        let starts_checksum = read_u32(bytes, &mut offset, "span starts checksum")?;
        let lengths_compressed_offset =
            read_u64(bytes, &mut offset, "span lengths compressed offset")?;
        let lengths_compressed_length =
            read_u32(bytes, &mut offset, "span lengths compressed length")?;
        let lengths_uncompressed_length =
            read_u32(bytes, &mut offset, "span lengths uncompressed length")?;
        let lengths_checksum = read_u32(bytes, &mut offset, "span lengths checksum")?;
        let start_encoding = read_u8(bytes, &mut offset, "span start encoding")?;
        let length_encoding = read_u8(bytes, &mut offset, "span length encoding")?;
        let starts_compression = read_u8(bytes, &mut offset, "span starts compression")?;
        let lengths_compression = read_u8(bytes, &mut offset, "span lengths compression")?;
        let reserved = read_u32(bytes, &mut offset, "span reserved")?;
        if span_count == 0
            || span_count > spans_per_block
            || reference_id >= reference_count
            || start_encoding != START_ENCODING_DELTA
            || length_encoding == 1
            || length_encoding > 2
            || starts_compression > 1
            || lengths_compression > 1
            || reserved != 0
        {
            return Err(Error::Corrupt("invalid span directory entry".into()));
        }
        let expected_first_id = if index == 0 { 0 } else { previous_span_id };
        if first_span_id != expected_first_id {
            return Err(Error::Corrupt(
                "span directory IDs are not contiguous".into(),
            ));
        }
        if let Some(previous_reference) = previous_reference
            && (reference_id < previous_reference
                || (reference_id == previous_reference && first_start < previous_start))
        {
            return Err(Error::Corrupt(
                "span directory is not coordinate ordered".into(),
            ));
        }
        let starts_end = starts_compressed_offset
            .checked_add(u64::from(starts_compressed_length))
            .ok_or_else(|| Error::Corrupt("span block range overflow".into()))?;
        if starts_compressed_offset < previous_starts_data_end || starts_end > starts_data_length {
            return Err(Error::Corrupt(
                "span starts block is outside data section".into(),
            ));
        }
        let lengths_end = lengths_compressed_offset
            .checked_add(u64::from(lengths_compressed_length))
            .ok_or_else(|| Error::Corrupt("span lengths range overflow".into()))?;
        if lengths_compressed_offset < previous_lengths_data_end
            || lengths_end > lengths_data_length
        {
            return Err(Error::Corrupt(
                "span lengths block is outside data section".into(),
            ));
        }
        if (span_count > 1 && starts_compressed_length == 0)
            || lengths_compressed_length == 0
            || lengths_uncompressed_length < LENGTH_PAYLOAD_HEADER_SIZE as u32
            || u64::from(starts_uncompressed_length) > MAX_BLOCK_BYTES
            || u64::from(lengths_uncompressed_length) > MAX_BLOCK_BYTES
        {
            return Err(Error::Corrupt("span block is too large".into()));
        }
        total_rows = total_rows
            .checked_add(u64::from(span_count))
            .ok_or_else(|| Error::Corrupt("span row count overflow".into()))?;
        previous_span_id = first_span_id
            .checked_add(u64::from(span_count))
            .ok_or_else(|| Error::Corrupt("span ID range overflow".into()))?;
        previous_starts_data_end = starts_end;
        previous_lengths_data_end = lengths_end;
        previous_reference = Some(reference_id);
        previous_start = first_start;
        entries.push(SpanDirectoryEntry {
            first_span_id,
            span_count,
            reference_id,
            first_start,
            starts_compressed_offset,
            starts_compressed_length,
            starts_uncompressed_length,
            starts_checksum,
            lengths_compressed_offset,
            lengths_compressed_length,
            lengths_uncompressed_length,
            lengths_checksum,
            start_encoding,
            length_encoding,
            starts_compression,
            lengths_compression,
        });
    }
    if total_rows != expected_span_count {
        return Err(Error::Corrupt("span row count mismatch".into()));
    }
    Ok(entries)
}

pub(crate) fn decode_for_stream(
    bytes: &[u8],
    offset: &mut usize,
    length: usize,
    count: usize,
    base: u64,
    bit_width: u8,
    context: &str,
) -> Result<Vec<u64>> {
    if bit_width > 63 {
        return Err(Error::Corrupt(format!("invalid {context} bit width")));
    }
    let bit_count = count
        .checked_mul(bit_width as usize)
        .ok_or_else(|| Error::Corrupt(format!("{context} bit count overflow")))?;
    let expected_length = bit_count
        .checked_add(7)
        .ok_or_else(|| Error::Corrupt(format!("{context} length overflow")))?
        / 8;
    if expected_length != length {
        return Err(Error::Corrupt(format!("{context} byte length mismatch")));
    }
    let end = offset
        .checked_add(length)
        .ok_or_else(|| Error::Corrupt(format!("{context} range overflow")))?;
    let packed = bytes
        .get(*offset..end)
        .ok_or_else(|| Error::Corrupt(format!("truncated {context}")))?;
    *offset = end;
    if bit_count % 8 != 0 {
        let used = bit_count % 8;
        let padding_mask = !((1_u8 << used) - 1);
        if packed.last().copied().unwrap_or(0) & padding_mask != 0 {
            return Err(Error::Corrupt(format!("nonzero padding bits in {context}")));
        }
    }
    let mut values = Vec::with_capacity(count);
    for index in 0..count {
        let bit_offset = index
            .checked_mul(bit_width as usize)
            .ok_or_else(|| Error::Corrupt(format!("{context} bit offset overflow")))?;
        let mut relative = 0_u64;
        for bit in 0..bit_width as usize {
            let target = bit_offset + bit;
            if packed[target / 8] & (1 << (target % 8)) != 0 {
                relative |= 1_u64 << bit;
            }
        }
        values.push(
            base.checked_add(relative)
                .ok_or_else(|| Error::Corrupt(format!("{context} value overflow")))?,
        );
    }
    Ok(values)
}

pub(crate) fn decode_delta_start_payload(
    payload: &[u8],
    entry: &SpanDirectoryEntry,
) -> Result<Vec<u64>> {
    let count = usize::try_from(entry.span_count)
        .map_err(|_| Error::Corrupt("span start count overflows usize".into()))?;
    if count == 0 || count > 100_000_000 {
        return Err(Error::Corrupt("invalid span start count".into()));
    }
    if payload.len() as u64 > MAX_BLOCK_BYTES {
        return Err(Error::Corrupt("span start payload is too large".into()));
    }
    let mut starts = Vec::with_capacity(count);
    starts.push(entry.first_start);
    let mut offset = 0_usize;
    let mut previous = entry.first_start;
    for _ in 1..count {
        let delta = read_canonical_varint(payload, &mut offset, "span start delta")?;
        previous = previous
            .checked_add(delta)
            .ok_or_else(|| Error::Corrupt("span start overflows coordinate".into()))?;
        starts.push(previous);
    }
    if offset != payload.len() {
        return Err(Error::Corrupt("trailing span start bytes".into()));
    }
    Ok(starts)
}

pub(crate) fn decode_length_payload(
    payload: &[u8],
    entry: &SpanDirectoryEntry,
) -> Result<Vec<u64>> {
    if payload.len() < LENGTH_PAYLOAD_HEADER_SIZE {
        return Err(Error::Corrupt("truncated span length payload".into()));
    }
    let mut offset = 0;
    let count = read_u32(payload, &mut offset, "span length count")?;
    let encoding = read_u8(payload, &mut offset, "span length encoding")?;
    let bit_width = read_u8(payload, &mut offset, "span length bit width")?;
    let reserved = read_u16(payload, &mut offset, "span length reserved")?;
    let base = read_u64(payload, &mut offset, "span length base")?;
    let stream_length = read_u32(payload, &mut offset, "span length stream length")? as usize;
    if reserved != 0
        || count == 0
        || count != entry.span_count
        || encoding != entry.length_encoding
        || encoding == 1
        || encoding > 2
    {
        return Err(Error::Corrupt("invalid span length metadata".into()));
    }
    if encoding == 0 && (bit_width != 0 || base != 0) {
        return Err(Error::Corrupt("invalid varint length metadata".into()));
    }
    if encoding == 2 && bit_width > 63 {
        return Err(Error::Corrupt("invalid FOR length bit width".into()));
    }
    let end = LENGTH_PAYLOAD_HEADER_SIZE
        .checked_add(stream_length)
        .ok_or_else(|| Error::Corrupt("span length range overflow".into()))?;
    if end != payload.len() {
        return Err(Error::Corrupt("span length stream length mismatch".into()));
    }
    let count = usize::try_from(count)
        .map_err(|_| Error::Corrupt("span length count overflows usize".into()))?;
    let stream = &payload[LENGTH_PAYLOAD_HEADER_SIZE..end];
    let lengths = match encoding {
        0 => {
            let mut stream_offset = 0;
            let mut values = Vec::with_capacity(count);
            for _ in 0..count {
                values.push(read_varint(stream, &mut stream_offset, "span length")?);
            }
            if stream_offset != stream.len() {
                return Err(Error::Corrupt("trailing span length bytes".into()));
            }
            values
        }
        2 => {
            let mut stream_offset = 0;
            let values = decode_for_stream(
                stream,
                &mut stream_offset,
                stream.len(),
                count,
                base,
                bit_width,
                "span length",
            )?;
            if stream_offset != stream.len() {
                return Err(Error::Corrupt("trailing span length bytes".into()));
            }
            values
        }
        _ => return Err(Error::Corrupt("unknown span length encoding".into())),
    };
    if lengths.contains(&0) {
        return Err(Error::Corrupt("span length is zero".into()));
    }
    Ok(lengths)
}

pub(crate) fn decode_span_rows(
    starts_payload: &[u8],
    lengths_payload: &[u8],
    entry: &SpanDirectoryEntry,
) -> Result<Vec<Span>> {
    let starts = decode_delta_start_payload(starts_payload, entry)?;
    let lengths = decode_length_payload(lengths_payload, entry)?;
    if starts.len() != lengths.len() {
        return Err(Error::Corrupt("span component row count mismatch".into()));
    }
    starts
        .into_iter()
        .zip(lengths)
        .map(|(start, length)| {
            start
                .checked_add(length)
                .ok_or(Error::InvalidCoordinate)
                .map(|_| Span {
                    reference_id: entry.reference_id,
                    start,
                    length,
                })
        })
        .collect()
}

#[cfg(test)]
pub(crate) fn decode_span_row(
    starts_payload: &[u8],
    lengths_payload: &[u8],
    entry: &SpanDirectoryEntry,
    row: usize,
) -> Result<Span> {
    decode_span_rows(starts_payload, lengths_payload, entry)?
        .get(row)
        .copied()
        .ok_or_else(|| Error::Corrupt("span row is out of bounds".into()))
}

/// An indexed GFF3 or BED source, its TBI/CSI coordinate index, and a GAI reader.
pub struct IndexedSource {
    source_path: PathBuf,
    source_format: SortFormat,
    coordinate_index_path: PathBuf,
    coordinate_index: CoordinateIndex,
    dictionary: CoordinateDictionary,
    reference_ids: HashMap<String, u32>,
    configured_terms: HashSet<String>,
    name_index: NameIndexReader,
}

impl IndexedSource {
    /// Opens a source, coordinate index, and GAI and rejects stale pairs by
    /// checking all source, index, and reference-dictionary fingerprints.
    pub fn open(
        source_path: impl AsRef<Path>,
        coordinate_index_path: impl AsRef<Path>,
        gai_path: impl AsRef<Path>,
    ) -> Result<Self> {
        Self::open_inner(
            source_path.as_ref(),
            coordinate_index_path.as_ref(),
            gai_path.as_ref(),
            None,
        )
    }

    /// Opens an indexed source while additionally checking the caller's
    /// attribute and normalization configuration against the stored header.
    pub fn open_with_options(
        source_path: impl AsRef<Path>,
        coordinate_index_path: impl AsRef<Path>,
        gai_path: impl AsRef<Path>,
        options: &NameIndexOptions,
    ) -> Result<Self> {
        let options = NameIndexOptions::new(options.attributes.clone(), options.case_sensitive)?;
        Self::open_inner(
            source_path.as_ref(),
            coordinate_index_path.as_ref(),
            gai_path.as_ref(),
            Some(&options),
        )
    }

    fn open_inner(
        source_path: &Path,
        coordinate_index_path: &Path,
        gai_path: &Path,
        options: Option<&NameIndexOptions>,
    ) -> Result<Self> {
        let source_path = source_path.to_path_buf();
        let coordinate_index_path = coordinate_index_path.to_path_buf();
        let source_format = SortFormat::from_path(&source_path)?;
        let name_index = NameIndexReader::open_mmap(gai_path)?;
        let (coordinate_index, coordinate_index_fingerprint) =
            read_coordinate_index_with_fingerprint(&coordinate_index_path)?;
        let dictionary = coordinate_index.dictionary(source_format)?;
        if fingerprint_file(&source_path)? != name_index.metadata.source_fingerprint {
            return Err(Error::Stale("source fingerprint does not match".into()));
        }
        if coordinate_index_fingerprint != name_index.metadata.coordinate_index_fingerprint {
            return Err(Error::Stale("TBI/CSI fingerprint does not match".into()));
        }
        if dictionary.fingerprint != name_index.metadata.reference_dictionary_fingerprint {
            return Err(Error::Stale(
                "reference dictionary fingerprint does not match".into(),
            ));
        }
        if name_index.metadata.reference_count
            != u32::try_from(dictionary.names.len())
                .map_err(|_| Error::Corrupt("reference dictionary is too large".into()))?
        {
            return Err(Error::Corrupt(
                "GAI reference count does not match coordinate index".into(),
            ));
        }
        if let Some(options) = options
            && (options.attributes != name_index.metadata.attributes
                || options.case_sensitive != name_index.metadata.case_sensitive)
        {
            return Err(Error::Stale(
                "configured attributes or normalization policy does not match GAI".into(),
            ));
        }
        for span_directory_entry in &name_index.span_directory {
            if usize::try_from(span_directory_entry.reference_id)
                .ok()
                .is_none_or(|reference_id| reference_id >= dictionary.names.len())
            {
                return Err(Error::Corrupt(
                    "span block has an invalid reference ID".into(),
                ));
            }
        }
        let reference_ids = dictionary
            .names
            .iter()
            .enumerate()
            .map(|(index, name)| (name.clone(), index as u32))
            .collect();
        let configured_terms = name_index.metadata.attributes.iter().cloned().collect();
        Ok(Self {
            source_path,
            source_format,
            coordinate_index_path,
            coordinate_index,
            dictionary,
            reference_ids,
            configured_terms,
            name_index,
        })
    }

    /// Returns GAI metadata.
    pub fn metadata(&self) -> &IndexMetadata {
        self.name_index.metadata()
    }

    /// Queries a configured attribute value or BED name exactly and returns matching
    /// records in source coordinate order. Each posted span is queried as an
    /// exact coordinate interval; returned chunks are then merged and read
    /// once. Unknown terms are successful empty queries.
    pub fn query_name(&mut self, term: &str) -> Result<Vec<GffRecord>> {
        self.query_name_with_mode(term, MatchMode::Exact)
    }

    /// Queries a configured attribute value exactly and returns records plus
    /// bounded query instrumentation.
    pub fn query_name_with_stats(&mut self, term: &str) -> Result<(Vec<GffRecord>, QueryStats)> {
        self.query_name_with_mode_and_stats(term, MatchMode::Exact)
    }

    /// Queries a configured GFF3 attribute value or BED name using an
    /// explicit match mode.
    /// Prefix mode streams matching FST keys through its prefix automaton;
    /// contains and regex modes scan keys incrementally. Matching spans are
    /// unioned and retrieved through the same batched TBI/CSI path. Regex
    /// patterns are unanchored searches by default, preserve regex syntax,
    /// and use Unicode-aware case folding for case-insensitive indexes.
    pub fn query_name_with_mode(
        &mut self,
        term: &str,
        match_mode: MatchMode,
    ) -> Result<Vec<GffRecord>> {
        self.query_name_with_mode_and_stats(term, match_mode)
            .map(|(records, _)| records)
    }

    /// Queries a configured attribute value or BED name using an explicit match mode and
    /// returns bounded query instrumentation.
    pub fn query_name_with_mode_and_stats(
        &mut self,
        term: &str,
        match_mode: MatchMode,
    ) -> Result<(Vec<GffRecord>, QueryStats)> {
        let matcher =
            CompiledMatch::new(term, match_mode, self.name_index.metadata.case_sensitive)?;
        if matcher.is_empty() {
            return Ok((Vec::new(), QueryStats::default()));
        }
        let span_ids = self
            .name_index
            .lookup_span_ids_with_matcher_and_stats(&matcher)?
            .0;
        let mut stats = QueryStats {
            requested_spans: span_ids.len() as u64,
            ..QueryStats::default()
        };
        if span_ids.is_empty() {
            return Ok((Vec::new(), stats));
        }

        let (spans, distinct_span_blocks, _) =
            self.name_index.resolve_span_ids_with_stats(&span_ids)?;
        stats.distinct_span_blocks_decoded = distinct_span_blocks;
        let requested_spans = spans
            .iter()
            .map(|span| SpanKey {
                reference_id: span.reference_id,
                start: span.start,
                length: span.length,
            })
            .collect::<HashSet<_>>();
        let mut raw_chunks = Vec::new();
        for span in &spans {
            let reference_sequence_id = usize::try_from(span.reference_id)
                .map_err(|_| Error::Corrupt("span reference ID overflows usize".into()))?;
            let reference_name = self
                .dictionary
                .names
                .get(reference_sequence_id)
                .ok_or_else(|| Error::Corrupt("span references invalid reference ID".into()))?;
            let query_start = Position::try_from(
                usize::try_from(span.start.checked_add(1).ok_or(Error::InvalidCoordinate)?)
                    .map_err(|_| Error::InvalidCoordinate)?,
            )
            .map_err(|_| Error::InvalidCoordinate)?;
            let query_end = Position::try_from(
                usize::try_from(span.end()?).map_err(|_| Error::InvalidCoordinate)?,
            )
            .map_err(|_| Error::InvalidCoordinate)?;
            let region = Region::new(reference_name.as_str(), query_start..=query_end);
            let chunks = match &self.coordinate_index {
                CoordinateIndex::Tabix(index) => {
                    index.query(reference_sequence_id, region.interval())?
                }
                CoordinateIndex::Csi(index) => {
                    index.query(reference_sequence_id, region.interval())?
                }
            };
            stats.exact_interval_queries += 1;
            stats.raw_chunks += chunks.len() as u64;
            raw_chunks.extend(chunks);
        }
        let merged_chunks = merge_query_chunks(raw_chunks);
        stats.merged_chunks = merged_chunks.len() as u64;
        let context = QueryReadContext {
            requested_spans: &requested_spans,
            reference_ids: &self.reference_ids,
            configured_terms: &self.configured_terms,
            case_sensitive: self.name_index.metadata.case_sensitive,
            matcher: &matcher,
        };
        let records = read_query_chunks(
            &self.source_path,
            self.source_format,
            &merged_chunks,
            &context,
            &mut stats,
        )?;
        Ok((records, stats))
    }

    /// Returns the path of the coordinate index used for this reader.
    pub fn coordinate_index_path(&self) -> &Path {
        &self.coordinate_index_path
    }

    /// Returns the GAI reader for callers needing direct span lookup.
    pub fn name_index(&self) -> &NameIndexReader {
        &self.name_index
    }
}

pub(crate) fn merge_query_chunks(mut chunks: Vec<Chunk>) -> Vec<Chunk> {
    chunks.sort_unstable_by_key(|chunk| (chunk.start(), chunk.end()));
    let mut merged: Vec<Chunk> = Vec::with_capacity(chunks.len());
    for chunk in chunks {
        if let Some(last) = merged.last_mut()
            && chunk.start() <= last.end()
        {
            if chunk.end() > last.end() {
                *last = Chunk::new(last.start(), chunk.end());
            }
            continue;
        }
        merged.push(chunk);
    }
    merged
}

pub(crate) struct QueryReadContext<'a> {
    pub(crate) requested_spans: &'a HashSet<SpanKey>,
    pub(crate) reference_ids: &'a HashMap<String, u32>,
    pub(crate) configured_terms: &'a HashSet<String>,
    pub(crate) case_sensitive: bool,
    pub(crate) matcher: &'a CompiledMatch,
}

pub(crate) fn read_query_chunks(
    source_path: &Path,
    source_format: SortFormat,
    chunks: &[Chunk],
    context: &QueryReadContext<'_>,
    stats: &mut QueryStats,
) -> Result<Vec<GffRecord>> {
    let source = File::open(source_path)?;
    match source_format {
        SortFormat::Gff => gff::read_gff_query_chunks(source, chunks, context, stats),
        SortFormat::Bed => bed::read_bed_query_chunks(source, chunks, context, stats),
    }
}
