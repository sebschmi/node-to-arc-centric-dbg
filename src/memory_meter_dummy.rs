use log::info;
use self_meter::Meter;
use std::time::Duration;

pub struct MemoryMeter;

impl MemoryMeter {
    pub fn new() -> Self {
        Self
    }

    pub fn report(&mut self) {
        info!("Memory reporting only supported on Linux");
    }
}
