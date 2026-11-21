//! Measure synchronous allocation audits without counting the libtest runner thread.
//! Take a ledger snapshot before the baseline so any lazy platform mutex setup
//! is included in fixture construction.

use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::Cell;

#[derive(Clone, Copy)]
struct Usage {
    current: isize,
    peak: isize,
}

thread_local! {
    // Constant initialization and no destructor keep this usable in allocation hooks.
    static USAGE: Cell<Usage> = const { Cell::new(Usage { current: 0, peak: 0 }) };
}

pub(crate) struct ThreadPeakAlloc;

impl ThreadPeakAlloc {
    pub(crate) fn current_usage(&self) -> isize {
        USAGE.with(|usage| usage.get().current)
    }

    pub(crate) fn peak_usage(&self) -> isize {
        USAGE.with(|usage| usage.get().peak)
    }

    pub(crate) fn reset_peak_usage(&self) {
        USAGE.with(|usage| {
            let current = usage.get().current;
            usage.set(Usage {
                current,
                peak: current,
            });
        });
    }
}

fn record(delta: isize) {
    let _ = USAGE.try_with(|usage| {
        let previous = usage.get();
        // Net bytes are signed because a thread can drop a value it didn't allocate.
        let current = previous.current.wrapping_add(delta);
        usage.set(Usage {
            current,
            peak: previous.peak.max(current),
        });
    });
}

// SAFETY: Allocation and deallocation delegate to System with the caller's
// pointer and layout. The counters do not allocate or unwind. Default realloc
// measures the old and new blocks together, as the production meter does.
unsafe impl GlobalAlloc for ThreadPeakAlloc {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        let ptr = unsafe { System.alloc(layout) };
        if !ptr.is_null() {
            record(layout.size() as isize);
        }
        ptr
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        unsafe { System.dealloc(ptr, layout) };
        record(-(layout.size() as isize));
    }
}
