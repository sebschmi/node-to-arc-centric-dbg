use log::info;
use self_meter::Meter;
use std::time::Duration;

pub struct MemoryMeter {
    meter: Meter,
}

impl MemoryMeter {
    pub fn new() -> Self {
        let mut meter = Meter::new(Duration::from_secs(1)).unwrap();
        meter.track_current_thread("main");
        meter.scan().unwrap();
        Self { meter }
    }

    pub fn report(&mut self) {
        self.meter.scan().unwrap();
        info!(
            "Current memory usage: {:.0}MiB",
            self.meter.report().unwrap().memory_rss as f64 / (1024.0 * 1024.0)
        );
    }
}
