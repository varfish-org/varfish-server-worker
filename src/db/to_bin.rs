//! Code for converting text-based files to binary ones (via flatbuffers).

/// Code for converting gene region from text-based to binary format.
pub mod gene_region {
    use std::{
        fs::File,
        io::Write,
        path::{Path, PathBuf},
        time::Instant,
    };

    use thousands::Separable;
    use tracing::debug;

    use crate::{
        common::{
            build_chrom_map, md5sum, numeric_gene_id, open_read_maybe_gz, sha256sum, trace_rss_now,
        },
        db::conf::DbDef,
        world_flatbuffers::var_fish_server_worker::{
            GeneRegionDatabase, GeneRegionDatabaseArgs, GeneRegionRecord,
        },
    };

    /// Module with code supporting the parsing.
    mod input {
        use serde::Deserialize;

        /// Record as created by VarFish DB Downloader.
        #[derive(Debug, Deserialize)]
        pub struct Record {
            /// Chromosome name
            pub chromosome: String,
            /// 1-based start position
            pub begin: u32,
            /// 1-based end position
            pub end: u32,
            /// ENSEMBL or Entrez gene ID
            pub gene_id: String,
        }
    }

    /// Perform conversion to flatbuffers `.bin` file.
    pub fn convert_to_bin<P, Q>(
        path_input_tsv: P,
        path_output_bin: Q,
        path_worker_db: &PathBuf,
        dbdef: &mut DbDef,
    ) -> Result<(), anyhow::Error>
    where
        P: AsRef<Path>,
        Q: AsRef<Path>,
    {
        debug!(
            "Converting gene regions from BED {:?} to binary {:?}",
            path_input_tsv.as_ref(),
            path_output_bin.as_ref()
        );
        let chrom_map = build_chrom_map();

        // Setup CSV reader for BED file - header is written as comment and must be ignored.
        let mut reader = csv::ReaderBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .comment(Some(b'#'))
            .from_reader(open_read_maybe_gz(path_input_tsv.as_ref())?);
        let before_parsing = Instant::now();

        let mut output_records = Vec::new();
        for record in reader.deserialize() {
            let record: input::Record = record?;
            output_records.push(GeneRegionRecord::new(
                *chrom_map
                    .get(&record.chromosome)
                    .unwrap_or_else(|| panic!("unknown chrom {:?}", &record.chromosome))
                    as u8,
                record.begin,
                record.end,
                numeric_gene_id(&record.gene_id)?,
            ));
        }

        debug!(
            "total time spent reading {:?} records: {:?}",
            output_records.len().separate_with_commas(),
            before_parsing.elapsed()
        );
        trace_rss_now();

        let before_writing = Instant::now();
        let mut builder = flatbuffers::FlatBufferBuilder::new();
        let records = builder.create_vector(output_records.as_slice());
        let gene_region_db = GeneRegionDatabase::create(
            &mut builder,
            &GeneRegionDatabaseArgs {
                records: Some(records),
            },
        );
        builder.finish_minimal(gene_region_db);
        let mut output_file = File::create(&path_output_bin)?;
        output_file.write_all(builder.finished_data())?;
        output_file.flush()?;
        debug!(
            "total time spent writing {} records: {:?}",
            output_records.len().separate_with_commas(),
            before_writing.elapsed()
        );

        // Update the DbDef.
        dbdef.bin_path = Some(
            path_output_bin
                .as_ref()
                .strip_prefix(path_worker_db)?
                .to_str()
                .unwrap()
                .to_string(),
        );
        dbdef.bin_md5 = Some(md5sum(path_output_bin.as_ref())?);
        dbdef.bin_sha256 = Some(sha256sum(path_output_bin.as_ref())?);

        Ok(())
    }
}

/// Code for converting xlink from text-based to binary format.
pub mod xlink {
    use std::{
        fs::File,
        io::Write,
        path::{Path, PathBuf},
        time::Instant,
    };

    use thousands::Separable;
    use tracing::debug;

    use crate::{
        common::{md5sum, numeric_gene_id, open_read_maybe_gz, sha256sum, trace_rss_now},
        db::conf::DbDef,
        world_flatbuffers::var_fish_server_worker::{
            XlinkDatabase, XlinkDatabaseArgs, XlinkRecord,
        },
    };

    /// Module with code for parsing the TSVs.
    pub mod input {
        use serde::Deserialize;

        #[derive(Debug, Deserialize)]
        pub struct Record {
            pub ensembl_gene_id: Option<String>,
            pub entrez_id: Option<u32>,
            pub gene_symbol: Option<String>,
        }
    }

    /// Perform conversion to flatbuffers `.bin` file.
    pub fn convert_to_bin<P, Q>(
        path_input_tsv: P,
        path_output_bin: Q,
        path_worker_db: &PathBuf,
        dbdef: &mut DbDef,
    ) -> Result<(), anyhow::Error>
    where
        P: AsRef<Path>,
        Q: AsRef<Path>,
    {
        let mut output_records = Vec::new();
        let mut output_strings = Vec::new();

        let before_parsing = Instant::now();

        let mut builder = flatbuffers::FlatBufferBuilder::new();

        debug!("parsing xlink TSV file from {:?}", path_input_tsv.as_ref());
        let mut reader = csv::ReaderBuilder::new()
            .has_headers(true)
            .delimiter(b'\t')
            .from_reader(open_read_maybe_gz(path_input_tsv)?);
        for record in reader.deserialize() {
            let record: input::Record = record?;
            if let (Some(entrez_id), Some(ensembl_gene_id), Some(gene_symbol)) =
                (record.entrez_id, record.ensembl_gene_id, record.gene_symbol)
            {
                output_records.push(XlinkRecord::new(
                    entrez_id,
                    numeric_gene_id(&ensembl_gene_id)?,
                ));
                let flat_str = builder.create_shared_string(&gene_symbol);
                output_strings.push(flat_str);
            }
        }

        debug!(
            "total time spent reading {} records: {:?}",
            output_records.len().separate_with_commas(),
            before_parsing.elapsed()
        );
        trace_rss_now();

        let before_writing = Instant::now();
        let records = builder.create_vector(&output_records);
        let strings = builder.create_vector(&output_strings);
        let xlink_db = XlinkDatabase::create(
            &mut builder,
            &XlinkDatabaseArgs {
                records: Some(records),
                symbols: Some(strings),
            },
        );
        builder.finish_minimal(xlink_db);
        let mut output_file = File::create(path_output_bin.as_ref())?;
        output_file.write_all(builder.finished_data())?;
        output_file.flush()?;
        debug!(
            "total time spent writing {} records: {:?}",
            output_records.len().separate_with_commas(),
            before_writing.elapsed()
        );
        trace_rss_now();

        // Update the DbDef.
        dbdef.bin_path = Some(
            path_output_bin
                .as_ref()
                .strip_prefix(path_worker_db)?
                .to_str()
                .unwrap()
                .to_string(),
        );
        dbdef.bin_md5 = Some(md5sum(path_output_bin.as_ref())?);
        dbdef.bin_sha256 = Some(sha256sum(path_output_bin.as_ref())?);

        Ok(())
    }
}

/// Code for converting ClinVar database to binary.
pub mod clinvar {
    use std::{
        fs::File,
        io::Write,
        path::{Path, PathBuf},
        time::Instant,
    };

    use thousands::Separable;
    use tracing::debug;

    use crate::{
        common::{build_chrom_map, md5sum, open_read_maybe_gz, sha256sum, trace_rss_now},
        db::conf::DbDef,
        world_flatbuffers::var_fish_server_worker::{
            ClinvarSvDatabase, ClinvarSvDatabaseArgs, ClinvarSvRecord,
        },
    };

    /// Module with code supporting the parsing.
    mod input {
        use std::collections::HashMap;

        use crate::sv::query::schema::{Pathogenicity, VariationType};
        use serde::{de, Deserialize, Deserializer};
        use tracing::warn;

        lazy_static::lazy_static! {
            static ref VARIATION_TYPE_LABELS: HashMap<&'static str, VariationType> = {
                let mut m = HashMap::new();
                m.insert("Complex", VariationType::Complex);
                m.insert("copy number gain", VariationType::Dup);
                m.insert("copy number loss", VariationType::Del);
                m.insert("Deletion", VariationType::Del);
                m.insert("Duplication", VariationType::Dup);
                m.insert("fusion", VariationType::Bnd);
                m.insert("Indel", VariationType::Cnv);
                m.insert("Insertion", VariationType::Ins);
                m.insert("Inversion", VariationType::Inv);
                m.insert("Microsatellite", VariationType::Microsatellite);
                m.insert("Tandem duplication", VariationType::Dup);
                m.insert("Translocation", VariationType::Bnd);
                m
            };
        }

        impl VariationType {
            pub fn from_label(label: &str) -> Result<VariationType, anyhow::Error> {
                if let Some(result) = VARIATION_TYPE_LABELS.get(label) {
                    Ok(*result)
                } else {
                    Err(anyhow::anyhow!("Invalid VariationType label: {}", label))
                }
            }
        }

        lazy_static::lazy_static! {
            static ref PATHOGENICITY_LABELS: HashMap<&'static str, Pathogenicity> = {
                let mut m = HashMap::new();
                m.insert("{\"benign\"}", Pathogenicity::Benign);
                m.insert("{\"benign\",\"likely benign\"}", Pathogenicity::LikelyBenign);
                m.insert("{\"likely benign\"}", Pathogenicity::LikelyBenign);
                m.insert("{\"likely pathogenic\"}", Pathogenicity::LikelyPathogenic);
                m.insert("{\"likely pathogenic\",\"pathogenic\"}", Pathogenicity::LikelyPathogenic);
                m.insert("{\"pathogenic\"}", Pathogenicity::Pathogenic);
                m.insert("{\"uncertain significance\"}", Pathogenicity::Uncertain);

                m
            };
        }

        impl Pathogenicity {
            pub fn from_label(label: &str) -> Result<Self, anyhow::Error> {
                if let Some(pathogenicity) = PATHOGENICITY_LABELS.get(label) {
                    Ok(*pathogenicity)
                } else {
                    warn!("Cannot decode pathogenicity from {}", label);
                    Ok(Pathogenicity::Uncertain)
                }
            }
        }

        /// Record as created by VarFish DB Downloader.
        #[derive(Debug, Deserialize)]
        pub struct Record {
            /// Chromosome name
            pub chromosome: String,
            /// 1-based start position
            pub begin: u32,
            /// 1-based end position
            pub end: u32,
            /// ClinVar SV variation type
            #[serde(deserialize_with = "from_variation_type_label")]
            pub variation_type: VariationType,
            /// The ClinVar VCV identifier
            pub vcv: String,
            /// Pathogenicity
            #[serde(
                alias = "summary_clinvar_pathogenicity",
                deserialize_with = "from_pathogenicity_summary"
            )]
            pub pathogenicity: Pathogenicity,
        }

        /// Deserialize "VariationType" from ClinVar TSV file
        ///
        /// This function will strip everything after the first underscore.
        fn from_variation_type_label<'de, D>(deserializer: D) -> Result<VariationType, D::Error>
        where
            D: Deserializer<'de>,
        {
            let s: &str = Deserialize::deserialize(deserializer)?;
            VariationType::from_label(s).map_err(de::Error::custom)
        }
        /// Deserialize "Pathogenicity" from ClinVar TSV file
        ///
        /// This function will strip everything after the first underscore.
        fn from_pathogenicity_summary<'de, D>(deserializer: D) -> Result<Pathogenicity, D::Error>
        where
            D: Deserializer<'de>,
        {
            let s: &str = Deserialize::deserialize(deserializer)?;
            Pathogenicity::from_label(s).map_err(de::Error::custom)
        }
    }

    /// Helper to convert VCV IDs to numbers.
    fn numeric_vcv_id(raw_id: &str) -> Result<u32, anyhow::Error> {
        let clean_id: String = raw_id
            .chars()
            .skip("VCV".len())
            .skip_while(|c| *c == '0')
            .collect();
        clean_id
            .parse::<u32>()
            .map_err(|e| anyhow::anyhow!("could not parse VCV id {:?}: {}", &clean_id, &e))
    }

    /// Perform conversion to flatbuffers `.bin` file.
    pub fn convert_to_bin<P, Q>(
        path_input_tsv: P,
        path_output_bin: Q,
        path_worker_db: &PathBuf,
        dbdef: &mut DbDef,
    ) -> Result<(), anyhow::Error>
    where
        P: AsRef<Path>,
        Q: AsRef<Path>,
    {
        let chrom_map = build_chrom_map();

        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(open_read_maybe_gz(path_input_tsv)?);
        let before_parsing = Instant::now();

        let mut output_records = Vec::new();
        for record in reader.deserialize() {
            let record: input::Record = record?;
            output_records.push(ClinvarSvRecord::new(
                *chrom_map.get(&record.chromosome).expect("unknown chrom") as u8,
                record.begin,
                record.end,
                record.variation_type.into(),
                record.pathogenicity.into(),
                numeric_vcv_id(&record.vcv)?,
            ));
        }

        debug!(
            "total time spent reading {} records: {:?}",
            output_records.len().separate_with_commas(),
            before_parsing.elapsed()
        );
        trace_rss_now();

        let before_writing = Instant::now();
        let mut builder = flatbuffers::FlatBufferBuilder::new();
        let records = builder.create_vector(output_records.as_slice());
        let clinvar_sv_db = ClinvarSvDatabase::create(
            &mut builder,
            &ClinvarSvDatabaseArgs {
                records: Some(records),
            },
        );
        builder.finish_minimal(clinvar_sv_db);
        let mut output_file = File::create(path_output_bin.as_ref())?;
        output_file.write_all(builder.finished_data())?;
        output_file.flush()?;
        debug!(
            "total time spent writing {} records: {:?}",
            output_records.len().separate_with_commas(),
            before_writing.elapsed()
        );

        // Update the DbDef.
        dbdef.bin_path = Some(
            path_output_bin
                .as_ref()
                .strip_prefix(path_worker_db)?
                .to_str()
                .unwrap()
                .to_string(),
        );
        dbdef.bin_md5 = Some(md5sum(path_output_bin.as_ref())?);
        dbdef.bin_sha256 = Some(sha256sum(path_output_bin.as_ref())?);

        Ok(())
    }
}

/// Code for converting other structural variant database to binary (incl. in-house).
pub mod vardbs {
    use std::fs::File;
    use std::io::Write;
    use std::path::{Path, PathBuf};
    use std::time::Instant;

    use anyhow::anyhow;
    use thousands::Separable;
    use tracing::debug;

    use crate::common::{build_chrom_map, md5sum, open_read_maybe_gz, sha256sum, trace_rss_now};
    use crate::db::conf::DbDef;
    use crate::db::inhouse::output::Record as InhouseDbRecord;
    use crate::sv::query::schema::SvType;
    use crate::world_flatbuffers::var_fish_server_worker::{
        BackgroundDatabase, BackgroundDatabaseArgs, BgDbRecord, SvType as FlatSvType,
    };

    use self::input::InputRecord;

    /// Code supporting the I/O of public database records and a common `InputRecord` for
    /// common representation.
    mod input {
        use serde::Deserialize;
        use tracing::error;

        use crate::db::inhouse::output::Record as InhouseDbRecord;
        use crate::sv::query::schema::SvType;

        /// dbVar database record as read from TSV file.
        #[derive(Debug, Deserialize)]
        pub struct DbVarRecord {
            /// chromosome name
            pub chromosome: String,
            /// begin position, 0-based
            pub begin: u32,
            /// end position, 0-based
            pub end: u32,
            /// number of overall carriers
            pub num_carriers: u32,
            /// type of the SV
            pub sv_type: String,
        }

        /// DGV database record as read from TSV file.
        #[derive(Debug, Deserialize)]
        pub struct DgvRecord {
            /// chromosome name
            pub chromosome: String,
            /// begin position, 0-based
            pub begin: u32,
            /// end position, 0-based
            pub end: u32,
            /// The structural variant type
            sv_type: String,
            /// Number of observed gains.
            observed_gains: u32,
            /// Number of observed losses
            observed_losses: u32,
        }

        /// DGV gold standard database record as read from TSV file.
        #[derive(Debug, Deserialize)]
        pub struct DgvGsRecord {
            /// chromosome name
            pub chromosome: String,
            /// begin position, 0-based
            pub begin_outer: u32,
            /// outer end position, 0-based
            pub end_outer: u32,
            /// The structural variant type
            pub sv_sub_type: String,
            /// Number of carriers.
            pub num_carriers: u32,
        }

        /// ExAC CNV database record as read from TSV file for deserialization from TSV.
        #[derive(Deserialize, Debug)]
        pub struct ExacRecord {
            /// chromosome name
            pub chromosome: String,
            /// outer start position, 0-based
            pub begin: u32,
            /// outer end position, 0-based
            pub end: u32,
            /// The structural vairant type
            pub sv_type: String,
        }

        /// Thousand Genomes SV database record as read from TSV file.
        #[derive(Debug, Deserialize)]
        pub struct G1kRecord {
            /// chromosome name
            pub chromosome: String,
            /// begin position, 0-based
            pub begin: u32,
            /// end position, 0-based
            pub end: u32,
            /// The structural vairant type
            pub sv_type: String,
            /// Number of hom. alt. alleles.
            pub n_homalt: u32,
            /// Number of het. alt. alleles.
            pub n_het: u32,
        }

        /// gnomAD SV database record as read from TSV file.
        #[derive(Debug, Deserialize)]
        pub struct GnomadRecord {
            /// chromosome name
            pub chromosome: String,
            /// begin position, 0-based
            pub begin: u32,
            /// end position, 0-based
            pub end: u32,
            /// The structural vairant type
            pub svtype: String,
            /// Number of homozygous alternative carriers
            pub n_homalt: u32,
            /// Number of heterozygous carriers
            pub n_het: u32,
        }

        /// Common type to convert input data to.
        pub struct InputRecord {
            /// Chromosome of start position.
            pub chromosome: String,
            /// Chromosome of end position.
            pub chromosome2: String,
            /// SV type
            pub sv_type: SvType,
            /// 0-based begin position
            pub begin: u32,
            /// 0-based end position
            pub end: u32,
            /// Number of carriers (or alleles), depending on database.
            pub count: u32,
        }

        impl TryInto<Option<InputRecord>> for InhouseDbRecord {
            type Error = &'static str;

            fn try_into(self) -> Result<Option<InputRecord>, Self::Error> {
                Ok(Some(InputRecord {
                    chromosome: self.chromosome,
                    chromosome2: self.chromosome2,
                    sv_type: self.sv_type,
                    begin: self.begin,
                    end: self.end,
                    count: self.carriers,
                }))
            }
        }

        impl TryInto<Option<InputRecord>> for DbVarRecord {
            type Error = &'static str;

            fn try_into(self) -> Result<Option<InputRecord>, Self::Error> {
                let sv_type = match self.sv_type.split(';').next().unwrap() {
                    "alu_insertion"
                    | "herv_insertion"
                    | "insertion"
                    | "line1_insertion"
                    | "mobile_element_insertion"
                    | "novel_sequence_insertion"
                    | "sva_insertion" => SvType::Ins,
                    "copy_number_gain" | "duplication" | "tandem_duplication" => SvType::Dup,
                    "alu_deletion" | "copy_number_loss" | "deletion" | "herv_deletion"
                    | "line1_deletion" | "sva_deletion" => SvType::Del,
                    "copy_number_variation" => SvType::Cnv,
                    _ => {
                        error!("sv_type = {}", &self.sv_type);
                        return Err("unknown SV type");
                    }
                };
                Ok(Some(InputRecord {
                    chromosome: self.chromosome.clone(),
                    chromosome2: self.chromosome,
                    begin: self.begin,
                    end: self.end,
                    sv_type,
                    count: 1,
                }))
            }
        }

        impl TryInto<Option<InputRecord>> for DgvRecord {
            type Error = &'static str;

            fn try_into(self) -> Result<Option<InputRecord>, Self::Error> {
                let sv_type = match self.sv_type.as_ref() {
                    "alu deletion"
                    | "deletion"
                    | "herv deletion"
                    | "line1 deletion"
                    | "mobile element deletion"
                    | "loss"
                    | "sva deletion" => SvType::Del,
                    "alu insertion"
                    | "herv insertion"
                    | "insertion"
                    | "line1 insertion"
                    | "mobile element insertion"
                    | "novel sequence insertion"
                    | "sva insertion" => SvType::Ins,
                    "duplication" | "gain" | "tandem duplication" => SvType::Dup,
                    "sequence alteration" | "complex" => return Ok(None), // skip
                    "gain+loss" | "CNV" => SvType::Cnv,
                    "inversion" => SvType::Inv,
                    "OTHER" => return Ok(None), // skip
                    _ => {
                        error!("sv_type = {}", &self.sv_type);
                        return Err("unknown SV type");
                    }
                };
                Ok(Some(InputRecord {
                    chromosome: self.chromosome.clone(),
                    chromosome2: self.chromosome,
                    begin: self.begin,
                    end: self.end,
                    sv_type,
                    count: self.observed_gains + self.observed_losses,
                }))
            }
        }

        impl TryInto<Option<InputRecord>> for DgvGsRecord {
            type Error = &'static str;

            fn try_into(self) -> Result<Option<InputRecord>, Self::Error> {
                let sv_type = match self.sv_sub_type.as_ref() {
                    "Gain" => SvType::Dup,
                    "Loss" => SvType::Del,
                    _ => {
                        error!("sv_type = {}", &self.sv_sub_type);
                        return Err("unknown SV type");
                    }
                };
                Ok(Some(InputRecord {
                    chromosome: self.chromosome.clone(),
                    chromosome2: self.chromosome,
                    begin: self.begin_outer,
                    end: self.end_outer,
                    sv_type,
                    count: self.num_carriers,
                }))
            }
        }

        impl TryInto<Option<InputRecord>> for ExacRecord {
            type Error = &'static str;

            fn try_into(self) -> Result<Option<InputRecord>, Self::Error> {
                let sv_type = match self.sv_type.as_ref() {
                    "DUP" => SvType::Dup,
                    "DEL" => SvType::Del,
                    _ => {
                        error!("sv_type = {}", &self.sv_type);
                        return Err("unknown SV type");
                    }
                };
                Ok(Some(InputRecord {
                    chromosome: self.chromosome.clone(),
                    chromosome2: self.chromosome,
                    begin: self.begin,
                    end: self.end,
                    sv_type,
                    count: 1,
                }))
            }
        }

        impl TryInto<Option<InputRecord>> for GnomadRecord {
            type Error = &'static str;

            fn try_into(self) -> Result<Option<InputRecord>, Self::Error> {
                let sv_type = match self.svtype.as_str() {
                    "CPX" => return Ok(None), // no correspondence
                    "CTX" | "BND" => SvType::Bnd,
                    "DEL" => SvType::Del,
                    "DUP" => SvType::Dup,
                    "INS" => SvType::Ins,
                    "INV" => SvType::Inv,
                    "MCNV" => SvType::Cnv,
                    _ => {
                        error!("sv_type = {}", &self.svtype);
                        return Err("unknown SV type");
                    }
                };
                Ok(Some(InputRecord {
                    chromosome: self.chromosome.clone(),
                    chromosome2: self.chromosome,
                    begin: self.begin,
                    end: self.end,
                    sv_type,
                    count: self.n_homalt + self.n_het,
                }))
            }
        }

        impl TryInto<Option<InputRecord>> for G1kRecord {
            type Error = &'static str;

            fn try_into(self) -> Result<Option<InputRecord>, Self::Error> {
                let sv_type = match self.sv_type.as_str() {
                    "CN0" | "CNV" => SvType::Cnv,
                    "DEL" => SvType::Del,
                    "DEL_ALU" | "DEL_HERV" | "DEL_LINE1" | "DEL_SVA" => SvType::Del,
                    "DUP" => SvType::Dup,
                    "INV" => SvType::Inv,
                    "INS" | "INS:ME:ALU" | "INS:ME:LINE1" | "INS:ME:SVA" => SvType::Ins,
                    _ => {
                        error!("sv_type = {}", &self.sv_type);
                        return Err("unknown SV type");
                    }
                };
                Ok(Some(InputRecord {
                    chromosome: self.chromosome.clone(),
                    chromosome2: self.chromosome,
                    begin: self.begin,
                    end: self.end,
                    sv_type,
                    count: self.n_homalt + self.n_het,
                }))
            }
        }
    }

    /// Known input file types.
    #[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
    pub enum InputFileType {
        Dbvar,
        Dgv,
        DgvGs,
        Exac,
        G1k,
        Gnomad,
        InhouseDb,
    }
    /// Deserialize from CSV reader to an `Option<records::InputRecord>`
    fn deserialize_loop<Rec>(
        reader: &mut csv::Reader<Box<dyn std::io::Read>>,
    ) -> Result<Vec<BgDbRecord>, anyhow::Error>
    where
        Rec: TryInto<Option<InputRecord>> + for<'de> serde::Deserialize<'de>,
        <Rec as TryInto<std::option::Option<InputRecord>>>::Error: core::fmt::Debug,
        <Rec as TryInto<std::option::Option<InputRecord>>>::Error: Send,
        <Rec as TryInto<std::option::Option<InputRecord>>>::Error: std::marker::Sync,
    {
        let chrom_map = build_chrom_map();
        let mut result = Vec::new();

        for record in reader.deserialize() {
            let record: Rec = record?;
            let maybe_record: Option<InputRecord> = record
                .try_into()
                .map_err(|err| anyhow!("problem with parsing: {:?}", &err))?;
            if let Some(record) = maybe_record {
                result.push(BgDbRecord::new(
                    *chrom_map
                        .get(&record.chromosome)
                        .unwrap_or_else(|| panic!("unknown chrom: {:?}", &record.chromosome))
                        as u8,
                    *chrom_map
                        .get(&record.chromosome2)
                        .unwrap_or_else(|| panic!("unknown chrom2: {:?}", &record.chromosome2))
                        as u8,
                    match record.sv_type {
                        SvType::Del => FlatSvType::Del,
                        SvType::Dup => FlatSvType::Dup,
                        SvType::Inv => FlatSvType::Inv,
                        SvType::Ins => FlatSvType::Ins,
                        SvType::Bnd => FlatSvType::Bnd,
                        SvType::Cnv => FlatSvType::Cnv,
                    },
                    record.begin,
                    record.end,
                    record.count,
                ));
            }
        }

        Ok(result)
    }

    /// Perform conversion to flatbuffers `.bin` file.
    pub fn convert_to_bin<P, Q>(
        path_input_tsv: P,
        path_output_bin: Q,
        path_worker_db: &PathBuf,
        input_type: InputFileType,
        dbdef: &mut DbDef,
    ) -> Result<(), anyhow::Error>
    where
        P: AsRef<Path>,
        Q: AsRef<Path>,
    {
        // Setup CSV reader for BED file - header is written as comment and must be ignored.
        let mut reader = csv::ReaderBuilder::new()
            .has_headers(false)
            .comment(Some(b'#'))
            .delimiter(b'\t')
            .from_reader(open_read_maybe_gz(path_input_tsv.as_ref())?);
        let before_parsing = Instant::now();

        let output_records = match input_type {
            InputFileType::Dbvar => deserialize_loop::<input::DbVarRecord>(&mut reader)?,
            InputFileType::Dgv => deserialize_loop::<input::DgvRecord>(&mut reader)?,
            InputFileType::DgvGs => deserialize_loop::<input::DgvGsRecord>(&mut reader)?,
            InputFileType::Exac => deserialize_loop::<input::ExacRecord>(&mut reader)?,
            InputFileType::G1k => deserialize_loop::<input::G1kRecord>(&mut reader)?,
            InputFileType::Gnomad => deserialize_loop::<input::GnomadRecord>(&mut reader)?,
            InputFileType::InhouseDb => deserialize_loop::<InhouseDbRecord>(&mut reader)?,
        };

        debug!(
            "total time spent reading {} records: {:?}",
            output_records.len().separate_with_commas(),
            before_parsing.elapsed()
        );
        trace_rss_now();

        let before_writing = Instant::now();
        let mut builder = flatbuffers::FlatBufferBuilder::new();
        let records = builder.create_vector(output_records.as_slice());
        let bg_db = BackgroundDatabase::create(
            &mut builder,
            &BackgroundDatabaseArgs {
                records: Some(records),
            },
        );
        builder.finish_minimal(bg_db);
        let mut output_file = File::create(path_output_bin.as_ref())?;
        output_file.write_all(builder.finished_data())?;
        output_file.flush()?;
        debug!(
            "total time spent writing {} records: {:?}",
            output_records.len().separate_with_commas(),
            before_writing.elapsed()
        );

        // Update the DbDef.
        dbdef.bin_path = Some(
            path_output_bin
                .as_ref()
                .strip_prefix(path_worker_db)?
                .to_str()
                .unwrap()
                .to_string(),
        );
        dbdef.bin_md5 = Some(md5sum(path_output_bin.as_ref())?);
        dbdef.bin_sha256 = Some(sha256sum(path_output_bin.as_ref())?);

        Ok(())
    }
}
