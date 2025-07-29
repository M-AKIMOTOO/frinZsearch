use std::fmt;
use std::io;
use std::error::Error;

#[derive(Debug)]
pub enum FrinZError {
    Io(io::Error),
    // Parse(String),
    Logic(String),
    Fit(String),
    // 必要に応じて他のエラータイプを追加
}

impl fmt::Display for FrinZError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            FrinZError::Io(err) => write!(f, "I/O Error: {}", err),
            // FrinZError::Parse(msg) => write!(f, "Parse Error: {}", msg),
            FrinZError::Logic(msg) => write!(f, "Logic Error: {}", msg),
            FrinZError::Fit(msg) => write!(f, "Fitting Error: {}", msg),
        }
    }
}

impl Error for FrinZError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            FrinZError::Io(err) => Some(err),
            _ => None,
        }
    }
}

impl From<io::Error> for FrinZError {
    fn from(err: io::Error) -> Self {
        FrinZError::Io(err)
    }
}

// 必要に応じて他のFrom実装を追加
// 例: std::num::ParseIntError, std::num::ParseFloatError など
