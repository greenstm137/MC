#include "mcPTBody.h"

mcPTBody::mcPTBody(void) {}
mcPTBody::~mcPTBody(void) {}

void mcPTBody::SetGrid(unsigned nx, unsigned ny, unsigned nz,
	double x0, double y0, double z0, double psx, double psy, double psz)
{
	nx_ = nx; ny_ = ny; nz_ = nz;
	size_ = nx * ny * nz;
	psx_ = psx; psy_ = psy; psz_ = psz;
	x0_ = x0; y0_ = y0; z0_ = z0;

	mediaIdxs_ = std::make_unique<std::vector<unsigned short>>(size_, 0);
	dens_ = std::make_unique<std::vector<double>>(size_, 0);

	// Hack!! Предполагаем что сред будет не больше 30.
	defaultDensities_ = std::make_unique<std::vector<double>>(30, 0);
}

void mcPTBody::SetDefaultDencities(const std::vector<double>& densities)
{
	defaultDensities_->resize(densities.size());
	for (auto i = 0; i < densities.size(); i++)
		(*defaultDensities_)[i] = densities[i];
}
