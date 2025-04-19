#include "mcScorePTGrid3D.h"
#include "mcGeometry.h"
#include "mcTransport.h"

mcScorePTGrid3D::mcScorePTGrid3D(const char* module_name, int nThreads, int nx, int ny, int nz)
	:mcScore(module_name, nThreads)
	, nx_(nx), ny_(nx), nz_(nz)
{
	energy_ = std::make_unique<std::vector<std::vector<double>>>(nThreads, std::vector<double>(nx * ny * nz, 0));
}

mcScorePTGrid3D::~mcScorePTGrid3D()
{
}

void mcScorePTGrid3D::ScoreFluence(const mcParticle& particle)
{
}

void mcScorePTGrid3D::ScorePoint(double edep
	, int iThread
	, const mcRegionReference& region
	, mc_particle_t pt
	, const geomVector3D& p)
{
	energy_->at(iThread)[region.idx_] += edep;
	etotal_[iThread] += edep;
}

void mcScorePTGrid3D::ScoreLine(double edep
	, int iThread
	, const mcRegionReference& region
	, mc_particle_t pt
	, const geomVector3D& p0
	, const geomVector3D& p1)
{
	energy_->at(iThread)[region.idx_] += edep;
	etotal_[iThread] += edep;
}

void mcScorePTGrid3D::CE2D()
{
	throw std::exception("mcScorePTGrid3D::CE2D: Not implemented");
}

void mcScorePTGrid3D::ClearEnergyMatrix()
{
	if (energy_ != nullptr)
	{
		for (int i = 0; i < nThreads_; i++)
		{
			memset(&(*energy_)[i][0], 0, nx_ * ny_ * nz_ * sizeof(double));
		}
	}
}

const std::vector<std::vector<double>>& mcScorePTGrid3D::EnergyMatrix() const
{
	return *energy_;
}

void mcScorePTGrid3D::dumpStatistic(ostream& os) const
{
}
