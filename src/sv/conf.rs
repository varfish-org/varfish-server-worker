//! Definition of records used by the ``sv`` sub command.

use serde::{Deserialize, Serialize};

/// Configuration for the database backing the SV annotation.
#[derive(Serialize, Deserialize, PartialEq, Debug)]
pub struct DbConf {
    /// Configuration of public databases.
    pub public_dbs: PublicDbsConf,
}

/// Configuration for public databases.
#[derive(Serialize, Deserialize, PartialEq, Debug)]
pub struct PublicDbsConf {
    /// Radius around BND sites used when building the database.
    pub slack_bnd: u32,
    /// Radius around INS sites used when building the database.
    pub slack_ins: u32,
    /// Minimal reciprocal overlap for SVs of the same type, used when building the database.
    pub min_overlap: f32,

    /// Relative path to the file with gnomAD-SV database with checksum.
    pub gnomad_sv: PathAndChecksum,
    /// Relative path to the file with dbVar SV database with checksum.
    pub dbvar: PathAndChecksum,
    /// Relative path to the file with DGV database with checksum.
    pub dgv: PathAndChecksum,
    /// Relative path to the file with DGV GS database with checksum.
    pub dgv_gs: PathAndChecksum,
    /// Relative path to the file with ExAC database with checksum.
    pub exac: PathAndChecksum,
}

/// A relative path to a file with its checksum.
#[derive(Serialize, Deserialize, PartialEq, Debug)]
pub struct PathAndChecksum {
    /// The (relative) path.
    pub path: String,
    /// Optional MD5 checksum.
    pub md5: Option<String>,
    /// Optional SHA256 checksum.
    pub sha256: Option<String>,
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;

    use crate::sv::conf::{PathAndChecksum, PublicDbsConf};

    use super::DbConf;

    #[test]
    fn test_parse_config_full() -> Result<(), anyhow::Error> {
        let toml_path = "tests/sv/conf/full.toml";
        let toml_str = std::fs::read_to_string(toml_path)?;
        let toml_data: DbConf = toml::from_str(&toml_str)?;

        assert_eq!(
            toml_data,
            DbConf {
                public_dbs: PublicDbsConf {
                    slack_bnd: 50,
                    slack_ins: 50,
                    min_overlap: 0.8,
                    gnomad_sv: PathAndChecksum {
                        path: "public/gnomad_sv.bin".to_owned(),
                        md5: Some("d41d8cd98f00b204e9800998ecf8427e".to_owned()),
                        sha256: Some(
                            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
                                .to_owned()
                        ),
                    },
                    dbvar: PathAndChecksum {
                        path: "public/dbvar.bin".to_owned(),
                        md5: Some("d41d8cd98f00b204e9800998ecf8427e".to_owned()),
                        sha256: Some(
                            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
                                .to_owned()
                        )
                    },
                    dgv: PathAndChecksum {
                        path: "public/dgv.bin".to_owned(),
                        md5: Some("d41d8cd98f00b204e9800998ecf8427e".to_owned()),
                        sha256: Some(
                            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
                                .to_owned()
                        ),
                    },
                    dgv_gs: PathAndChecksum {
                        path: "public/dgv_gs.bin".to_owned(),
                        md5: Some("d41d8cd98f00b204e9800998ecf8427e".to_owned()),
                        sha256: Some(
                            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
                                .to_owned()
                        ),
                    },
                    exac: PathAndChecksum {
                        path: "public/exac.bin".to_owned(),
                        md5: Some("d41d8cd98f00b204e9800998ecf8427e".to_owned()),
                        sha256: Some(
                            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
                                .to_owned()
                        ),
                    }
                }
            }
        );

        Ok(())
    }
}
