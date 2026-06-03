#include "openmc/photon_track.h"

#include <algorithm>
#include <string>

#include <fmt/format.h>

#include "openmc/bank.h"
#include "openmc/bank_io.h"
#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/hdf5_interface.h"
#include "openmc/material.h"
#include "openmc/mcpl_interface.h"
#include "openmc/message_passing.h"
#include "openmc/nuclide.h"
#include "openmc/output.h"
#include "openmc/particle.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/universe.h"

#ifdef OPENMC_MPI
#include <mpi.h>
#endif

namespace openmc {

namespace {

// hid_t h5_collision_track_banktype()
// {
//   hid_t banktype = H5Tcreate(H5T_COMPOUND, sizeof(CollisionTrackSite));
//   H5Tinsert(
//     banktype, "time", HOFFSET(CollisionTrackSite, time), H5T_NATIVE_DOUBLE);
//   H5Tinsert(banktype, "event_mt", HOFFSET(CollisionTrackSite, event_mt),
//     H5T_NATIVE_INT);
//   H5Tinsert(
//     banktype, "cell_id", HOFFSET(CollisionTrackSite, cell_id), H5T_NATIVE_INT);
//   H5Tinsert(banktype, "particle", HOFFSET(CollisionTrackSite, particle),
//     H5T_NATIVE_INT);
//   H5Tinsert(banktype, "parent_id", HOFFSET(CollisionTrackSite, parent_id),
//     H5T_NATIVE_INT64);
//   return banktype;
// }

hid_t h5_photon_track_banktype()
{
  hid_t banktype = H5Tcreate(H5T_COMPOUND, sizeof(PhotonTrackSite));
  H5Tinsert(banktype, "batch_no", HOFFSET(PhotonTrackSite, batch_no),H5T_NATIVE_INT64);
  H5Tinsert(banktype, "parent_id", HOFFSET(PhotonTrackSite, parent_id),H5T_NATIVE_INT64);
  H5Tinsert(banktype, "x", HOFFSET(PhotonTrackSite, r.x), H5T_NATIVE_DOUBLE);
  H5Tinsert(banktype, "y", HOFFSET(PhotonTrackSite, r.y), H5T_NATIVE_DOUBLE);
  H5Tinsert(banktype, "z", HOFFSET(PhotonTrackSite, r.z), H5T_NATIVE_DOUBLE);
  H5Tinsert(banktype, "time", HOFFSET(PhotonTrackSite, time), H5T_NATIVE_DOUBLE);
  H5Tinsert(banktype, "dE", HOFFSET(PhotonTrackSite, dE), H5T_NATIVE_DOUBLE);
  H5Tinsert(banktype, "event_mt", HOFFSET(PhotonTrackSite, event_mt), H5T_NATIVE_INT);
  H5Tinsert(banktype, "cell_id", HOFFSET(PhotonTrackSite, cell_id), H5T_NATIVE_INT);

  return banktype;
}

void write_photon_track_bank(hid_t group_id,
  openmc::span<PhotonTrackSite> photon_track_bank,
  const openmc::vector<int64_t>& bank_index)
{
  hid_t banktype = h5_photon_track_banktype();
#ifdef OPENMC_MPI
  write_bank_dataset("photon_track_bank", group_id, photon_track_bank,
    bank_index, banktype, banktype, mpi::photon_track_site);
#else
  write_bank_dataset("photon_track_bank", group_id, photon_track_bank,
    bank_index, banktype, banktype);
#endif

  H5Tclose(banktype);
}

void write_h5_photon_track(const char* filename,
  openmc::span<PhotonTrackSite> photon_track_bank,
  const openmc::vector<int64_t>& bank_index)
{
#ifdef PHDF5
  bool parallel = true;
#else
  bool parallel = false;
#endif

  if (!filename)
    fatal_error("write_h5_photon_track filename needs a nonempty name.");

  std::string filename_(filename);
  const auto extension = get_file_extension(filename_);
  if (extension.empty()) {
    filename_.append(".h5");
  } else if (extension != "h5") {
    warning("write_h5_photon_track was passed a file extension differing "
            "from .h5, but an hdf5 file will be written.");
  }

  hid_t file_id;
  if (mpi::master || parallel) {
    file_id = file_open(filename_.c_str(), 'w', true);

    // Write filetype and version info
    write_attribute(file_id, "filetype", "photon_track");
    write_attribute(file_id, "version", VERSION_PHOTON_TRACK);
  }

  write_photon_track_bank(file_id, photon_track_bank, bank_index);

  if (mpi::master || parallel)
    file_close(file_id);
}

} // namespace

// bool should_record_event(int id_cell, int mt_event, const std::string& nuclide,
//   int id_universe, int id_material, double energy_loss, std::string particle_type)
// {
//   std::cout << id_cell << " | " << mt_event << " | " << nuclide  << " | " << id_universe << " | " << id_material << " | " << energy_loss << " | " << particle_type << std::endl;
//   auto matches_filter = [](const auto& filter_set, const auto& value) {
//     return filter_set.empty() || filter_set.count(value) > 0;
//   };

//   const auto& cfg = settings::photon_track_config;
//   // std::cout << "Particle type: " << particle_type.str() << " should be " << cfg.particle_types[0] << std::endl;
//   // ParticleType checkType = ParticleType(cfg.particle_types);
//   return simulation::current_batch > settings::n_inactive &&
//          !simulation::photon_track_bank.full() &&
//          matches_filter(cfg.cell_ids, id_cell) &&
//          matches_filter(cfg.mt_numbers, mt_event) &&
//          matches_filter(cfg.universe_ids, id_universe) &&
//          matches_filter(cfg.material_ids, id_material) &&
//          matches_filter(cfg.nuclides, nuclide) &&
//          (cfg.deposited_energy_threshold == 0 ||
//            cfg.deposited_energy_threshold < energy_loss) &&
//          matches_filter(cfg.particle_types, particle_type);
// }

void photon_track_reserve_bank()
{
  simulation::photon_track_bank.reserve(
    settings::photon_track_config.max_collisions);
}

void photon_track_flush_bank()
{
  const auto& cfg = settings::photon_track_config;
  if (simulation::pt_current_file > cfg.max_files)
    return;

  bool last_batch = (simulation::current_batch == settings::n_batches);
  if (!simulation::photon_track_bank.full() && !last_batch)
    return;

  auto size = simulation::photon_track_bank.size();
  if (size == 0 && !last_batch)
    return;

  auto photon_track_work_index = mpi::calculate_parallel_index_vector(size);
  openmc::span<PhotonTrackSite> photontrackbankspan(
    simulation::photon_track_bank.begin(), size);

  std::string ext = "h5";
  auto filename = fmt::format("{}photon_track.{}.{}", settings::path_output,
    simulation::pt_current_file, ext);

  if (cfg.max_files == 1 || (simulation::pt_current_file == 1 && last_batch)) {
    filename = settings::path_output + "photon_track." + ext;
  }
  write_message("Creating {}...", filename, 4);

  // Remove option to write as mcpl
  write_h5_photon_track(filename.c_str(), photontrackbankspan, photon_track_work_index);


  simulation::photon_track_bank.clear();
  if (!last_batch && cfg.max_files >= 1) {
    photon_track_reserve_bank();
  }
  ++simulation::pt_current_file;
}

void photon_track_record(Particle& particle)
{
  int cell_index = particle.lowest_coord().cell();
  int cell_id = 0;
  if (cell_index != C_NONE) {
    cell_id = model::cells[cell_index]->id_;

  }
  double delta_E = particle.E_last();
  PhotonTrackSite site;
  site.r = particle.r();
  site.dE = delta_E;
  site.time = particle.time();
  site.event_mt = particle.event_mt();
  site.cell_id = cell_id;
  site.parent_id = particle.id();
  site.batch_no = simulation::current_batch;
  simulation::photon_track_bank.thread_safe_append(site);
}

} // namespace openmc
