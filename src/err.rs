use std::{
    num::ParseIntError,
    process::{ExitCode, Termination},
};

#[derive(thiserror::Error, Debug, Clone)]
pub enum AppError {
    // #[error("Internal error")]
    // Internal,
    #[error("Invalid path")]
    InvalidPath,
}

impl Termination for AppError {
    fn report(self) -> ExitCode {
        match self {
            // AppError::Internal => ExitCode::from(1),
            AppError::InvalidPath => ExitCode::from(1),
        }
    }
}

#[derive(thiserror::Error, Debug, Clone)]
pub enum ArgError {
    #[error("Invalid format in interval")]
    IntervalInvalidFormat,
    #[error("Invalid integer coordinates in interval")]
    IntervalInvalidInts(#[from] ParseIntError),
}
