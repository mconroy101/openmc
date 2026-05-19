#ifndef OPENMC_PHOTON_TRACK_H
#define OPENMC_PHOTON_TRACK_H

#include <string>

namespace openmc {

class Particle;

//! Reserve space in the collision track bank according to user settings.
void photon_track_reserve_bank();

//! Write collision track data to disk when the bank is full or the batch ends.
void photon_track_flush_bank();

//! Record the current particle as a collision-track entry when applicable.
//!
//! \param particle Particle whose collision should be recorded if eligible
void photon_track_record(Particle& particle);

} // namespace openmc

#endif // OPENMC_PHOTON_TRACK_H
