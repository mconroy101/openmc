#include "openmc/tallies/filter_photon_origin.h"

#include "openmc/cell.h"

namespace openmc {

void PhotonOriginFilter::get_all_bins(
  const Particle& p, TallyEstimator estimator, FilterMatch& match) const
{
  auto search = map_.find(p.photon_origin_cell()); // CURRENTLY THIS DOES NOT MATCH BECAUSE p.photon_origin_cell() STORES THE CELL ID, WHEREAS map_ STORES THE CELL INDICES
  if (search != map_.end()) {
    match.bins_.push_back(search->second);
    match.weights_.push_back(1.0);
  }
}

std::string PhotonOriginFilter::text_label(int bin) const
{
  return "Photon Origin Cell " + std::to_string(model::cells[cells_[bin]]->id_);
}

} // namespace openmc
