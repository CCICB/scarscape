//! Statistical engines
//!
//! Stats engines consume normalized event streams and reference metadata
//! to compute summary statistics.
//!
//! They assume inputs are already validated and normalized.

use crate::error::{Error, Result};
use crate::model::SmallVariantCounts;
use seqlib::mutations::DnaSmallMutation;

/// Summarise Variant Type Counts
/// From an iterable stream of DnaSmallMutation counts
pub fn tally_small_variant_types<I>(variant_stream: I) -> SmallVariantCounts
where
    I: IntoIterator<Item = Result<DnaSmallMutation>>,
{
    let mut total: u64 = 0;
    let mut skipped: u64 = 0;

    for variant_result in variant_stream {
        total += 1;
        if variant_result.is_err() {
            skipped += 1;
        }
    }

    SmallVariantCounts { total }
}
