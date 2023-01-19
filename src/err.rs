use std::{
    num::ParseIntError,
    process::{ExitCode, Termination},
};

#[derive(thiserror::Error, Debug, Clone)]
pub enum AppError {
    #[allow(dead_code)]
    #[error("Unknown contig")]
    UnknownContig,
}

impl Termination for AppError {
    fn report(self) -> ExitCode {
        match self {
            AppError::UnknownContig => ExitCode::from(1),
        }
    }
}

#[derive(thiserror::Error, Debug, Clone)]
pub enum ArgError {
    #[allow(dead_code)]
    #[error("Invalid format in interval")]
    IntervalInvalidFormat,
    #[error("Invalid integer coordinates in interval")]
    IntervalInvalidInts(#[from] ParseIntError),
}
