//! Module with datastructures stored in RocksDB.

use byteorder::{ByteOrder, LittleEndian};

/// Genotype counts.
#[derive(Debug, Default, Clone)]
pub struct Counts {
    /// Number of hom. ref. carriers.
    pub count_homref: u32,
    /// Number of hemi. ref. carriers.
    pub count_hemiref: u32,
    /// Number of het. carriers.
    pub count_het: u32,
    /// Number of hom. alt. carriers.
    pub count_homalt: u32,
    /// Number of hemi. alt. carriers.
    pub count_hemialt: u32,
}

impl Counts {
    /// Convert to a byte vector.
    pub fn to_vec(&self) -> Vec<u8> {
        let mut buf = Vec::with_capacity(20);
        buf.extend_from_slice(&self.count_homref.to_le_bytes());
        buf.extend_from_slice(&self.count_hemiref.to_le_bytes());
        buf.extend_from_slice(&self.count_het.to_le_bytes());
        buf.extend_from_slice(&self.count_homalt.to_le_bytes());
        buf.extend_from_slice(&self.count_hemialt.to_le_bytes());
        buf
    }

    /// Convert from a byte vector.
    pub fn from_vec(buf: &[u8]) -> Self {
        Self {
            count_homref: LittleEndian::read_u32(&buf[0..4]),
            count_hemiref: LittleEndian::read_u32(&buf[4..8]),
            count_het: LittleEndian::read_u32(&buf[8..12]),
            count_homalt: LittleEndian::read_u32(&buf[12..16]),
            count_hemialt: LittleEndian::read_u32(&buf[16..20]),
        }
    }

    /// Aggregate other into self.
    pub fn aggregate(&mut self, other: Self) {
        self.count_homref += other.count_homref;
        self.count_hemiref += other.count_hemiref;
        self.count_het += other.count_het;
        self.count_homalt += other.count_homalt;
        self.count_hemialt += other.count_hemialt;
    }
}

/// Genotype in a carrier.
#[derive(Debug, Default, Clone, Copy, PartialOrd, Ord, PartialEq, Eq)]
pub enum Genotype {
    #[default]
    HomRef,
    HemiRef,
    Het,
    HomAlt,
    HemiAlt,
}

impl Genotype {
    /// Convert to a byte.
    pub fn to_byte(self) -> u8 {
        match self {
            Genotype::HomRef => 0,
            Genotype::HemiRef => 1,
            Genotype::Het => 2,
            Genotype::HomAlt => 3,
            Genotype::HemiAlt => 4,
        }
    }
}

/// Error for `TryFrom<u8>` for `Genotype`.
#[derive(Debug, Clone, thiserror::Error)]
pub enum GenotypeTryFromByteError {
    #[error("invalid genotype byte: {0}")]
    InvalidByte(u8),
}

impl TryFrom<u8> for Genotype {
    type Error = GenotypeTryFromByteError;

    fn try_from(byte: u8) -> Result<Self, Self::Error> {
        match byte {
            0 => Ok(Genotype::HomRef),
            1 => Ok(Genotype::HemiRef),
            2 => Ok(Genotype::Het),
            3 => Ok(Genotype::HomAlt),
            4 => Ok(Genotype::HemiAlt),
            _ => Err(GenotypeTryFromByteError::InvalidByte(byte)),
        }
    }
}

/// Store one carrier by UUID and index in the pedigree.
#[derive(Debug, Default, Clone, PartialOrd, Ord, PartialEq, Eq)]
pub struct Carrier {
    /// Case UUID.
    pub uuid: uuid::Uuid,
    /// Index of the carrier in the pedigree.
    pub index: u8,
    /// Genotype of the carrier.
    pub genotype: Genotype,
}

/// Carrier UUIDs.
///
/// We store the UUIDs serialized as byte vectors and the index of the carrier in the pedigree.
#[derive(Debug, Default, Clone)]
pub struct CarrierList {
    /// List of carrier UUIDs.
    pub carriers: Vec<Carrier>,
}

impl CarrierList {
    /// Convert to a byte vector.
    pub fn to_vec(&self) -> Vec<u8> {
        let mut buf = Vec::with_capacity(2 + 18 * self.carriers.len());
        buf.extend_from_slice(&(self.carriers.len() as u16).to_le_bytes());
        for carrier in &self.carriers {
            buf.extend_from_slice(&carrier.uuid.as_u128().to_le_bytes());
            buf.push(carrier.index);
            buf.push(carrier.genotype.to_byte());
        }
        buf
    }

    /// Order the carriers by UUID.
    pub fn sort(&mut self) {
        self.carriers.sort();
    }

    /// Aggregate other into self.
    pub fn aggregate(&mut self, mut other: Self) {
        self.carriers.append(&mut other.carriers);
        self.carriers.sort();
        self.carriers.dedup();
    }
}

impl TryFrom<&[u8]> for CarrierList {
    type Error = GenotypeTryFromByteError;

    fn try_from(buf: &[u8]) -> Result<Self, Self::Error> {
        let mut carriers = Vec::with_capacity((buf.len() - 2) / 18);
        let num_carriers = LittleEndian::read_u16(&buf[0..2]) as usize;
        for i in 0..num_carriers {
            let uuid =
                uuid::Uuid::from_u128(LittleEndian::read_u128(&buf[2 + 18 * i..2 + 18 * i + 16]));
            let index = buf[2 + 18 * i + 16];
            let genotype = Genotype::try_from(buf[2 + 18 * i + 17])?;
            carriers.push(Carrier {
                uuid,
                index,
                genotype,
            });
        }
        Ok(Self { carriers })
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_counts() {
        let counts = Counts {
            count_homref: 1,
            count_hemiref: 2,
            count_het: 3,
            count_homalt: 4,
            count_hemialt: 5,
        };

        let buf = counts.to_vec();
        insta::assert_debug_snapshot!(&buf);
        assert_eq!(buf.len(), 20);

        let counts2 = Counts::from_vec(&buf);
        insta::assert_debug_snapshot!(&counts2);
    }

    #[test]
    fn test_carrier_list() -> Result<(), anyhow::Error> {
        let carrier_list = CarrierList {
            carriers: vec![
                Carrier {
                    uuid: uuid::Uuid::parse_str("00000000-0000-0000-0000-000000000000").unwrap(),
                    index: 0,
                    genotype: Genotype::HomRef,
                },
                Carrier {
                    uuid: uuid::Uuid::parse_str("00000000-0000-0000-0000-000000000001").unwrap(),
                    index: 1,
                    genotype: Genotype::HemiAlt,
                },
            ],
        };

        let buf = carrier_list.to_vec();
        insta::assert_debug_snapshot!(&buf);
        assert_eq!(buf.len(), 38);

        let carrier_list2 = CarrierList::try_from(buf.as_slice())?;
        insta::assert_debug_snapshot!(&carrier_list2);

        Ok(())
    }
}
