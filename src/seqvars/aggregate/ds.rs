//! Module with datastructures stored in RocksDB.

use byteorder::{ByteOrder, LittleEndian};

/// Genotype counts.
#[derive(Debug, Default, Clone)]
pub struct Counts {
    /// Number of alleles in samples.
    pub count_an: u32,
    /// Number of het. carriers.
    pub count_het: u32,
    /// Number of hom. carriers.
    pub count_hom: u32,
    /// Number of hemi. carriers.
    pub count_hemi: u32,
}

impl Counts {
    /// Convert to a byte vector.
    pub fn to_vec(&self) -> Vec<u8> {
        let mut buf = Vec::with_capacity(16);
        buf.extend_from_slice(&self.count_an.to_le_bytes());
        buf.extend_from_slice(&self.count_het.to_le_bytes());
        buf.extend_from_slice(&self.count_hom.to_le_bytes());
        buf.extend_from_slice(&self.count_hemi.to_le_bytes());
        buf
    }

    /// Convert from a byte vector.
    pub fn from_vec(buf: &[u8]) -> Self {
        Self {
            count_an: LittleEndian::read_u32(&buf[0..4]),
            count_het: LittleEndian::read_u32(&buf[4..8]),
            count_hom: LittleEndian::read_u32(&buf[8..12]),
            count_hemi: LittleEndian::read_u32(&buf[12..16]),
        }
    }

    /// Aggregate other into self.
    pub fn aggregate(&mut self, other: Self) {
        self.count_an += other.count_an;
        self.count_het += other.count_het;
        self.count_hom += other.count_hom;
        self.count_hemi += other.count_hemi;
    }
}

/// Genotype in a carrier.
#[derive(Debug, Default, Clone, Copy, PartialOrd, Ord, PartialEq, Eq)]
pub enum Genotype {
    #[default]
    HomRef, // same as hemizygous reference
    Het,
    HomAlt,
    HemiAlt,
}

impl Genotype {
    /// Convert to a byte.
    pub fn to_byte(self) -> u8 {
        match self {
            Genotype::HomRef => 0,
            Genotype::Het => 1,
            Genotype::HomAlt => 2,
            Genotype::HemiAlt => 3,
        }
    }

    /// Convert from a byte.
    pub fn from_byte(byte: u8) -> Self {
        match byte {
            0 => Genotype::HomRef,
            1 => Genotype::Het,
            2 => Genotype::HomAlt,
            3 => Genotype::HemiAlt,
            _ => panic!("invalid genotype byte"),
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

    /// Convert from a byte vector.
    pub fn from_vec(buf: &[u8]) -> Self {
        let mut carriers = Vec::with_capacity((buf.len() - 2) / 18);
        let num_carriers = LittleEndian::read_u16(&buf[0..2]) as usize;
        for i in 0..num_carriers {
            let uuid =
                uuid::Uuid::from_u128(LittleEndian::read_u128(&buf[2 + 18 * i..2 + 18 * i + 16]));
            let index = buf[2 + 18 * i + 16];
            let genotype = Genotype::from_byte(buf[2 + 18 * i + 17]);
            carriers.push(Carrier {
                uuid,
                index,
                genotype,
            });
        }
        Self { carriers }
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

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn counts() {
        let counts = Counts {
            count_an: 1,
            count_het: 2,
            count_hom: 3,
            count_hemi: 4,
        };

        let buf = counts.to_vec();
        insta::assert_debug_snapshot!(&buf);
        assert_eq!(buf.len(), 16);

        let counts2 = Counts::from_vec(&buf);
        insta::assert_debug_snapshot!(&counts2);
    }

    #[test]
    fn carrier_list() {
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

        let carrier_list2 = CarrierList::from_vec(&buf);
        insta::assert_debug_snapshot!(&carrier_list2);
    }
}
