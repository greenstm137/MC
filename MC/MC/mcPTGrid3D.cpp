#include "mcPTGrid3D.h"
#include "mcGeometry.h"
#include "mcPTBody.h"
#include "mcMedia.h"
#include "mcPhysics.h"
#include "mcThread.h"
#include "mcScorePTGrid3D.h"
#include <float.h>

#define GEOM_EPSILON	1E-6

mcPTGrid3D::mcPTGrid3D(const geomVector3D& orgn, const geomVector3D& vz, const geomVector3D& vx, const mcPTBody& mcbody, int nThreads)
	: mcTransportPrism(orgn, vz, vx), mcbody_(mcbody), nThreads_(nThreads)
{
	setColor(0., 1., 0., 0.5);

	nx_ = mcbody_.nx_;
	ny_ = mcbody_.ny_;
	nz_ = mcbody_.nz_;

	psx_ = mcbody_.psx_;
	psy_ = mcbody_.psy_;
	psz_ = mcbody_.psz_;

	x0_ = mcbody_.x0_;
	y0_ = mcbody_.y0_;
	z0_ = mcbody_.z0_;

	size_ = nx_ * ny_ * nz_;
	setGeometry(psx_ * nx_, psy_ * ny_, psz_ * nz_, x0_, y0_, z0_);

	score_ = new mcScorePTGrid3D(this->getName(), nThreads, nx_, ny_, nz_);
}

mcPTGrid3D::~mcPTGrid3D(void)
{
	delete score_;
}

int mcPTGrid3D::getIdxAtPoint(const geomVector3D& p, short* pgidx) const
{
	double dx = p.x() - x0_, dy = p.y() - y0_, dz = p.z() - z0_;
	if (dx < 0 || dy < 0 || dz < 0 || dx >= ax_ || dy >= ay_ || dz >= az_)
		return -1;
	else
	{
		int i = int(dx / psx_), j = int(dy / psy_), k = int(dz / psz_);
		pgidx[0] = i; pgidx[1] = j; pgidx[2] = k;
		return (k * ny_ + j) * nx_ + i;
	}
}

double mcPTGrid3D::getDistanceInsideVoxel(const mcParticle& particle, short* gidxNext, int& idx)
{
	// Расстояния до граней параллелепипеда
	idx = 0;
	gidxNext[0] = particle.region.gidx_[0];
	gidxNext[1] = particle.region.gidx_[1];
	gidxNext[2] = particle.region.gidx_[2];

	double dx = DBL_MAX, dy = DBL_MAX, dz = DBL_MAX;
	double vx = particle.u.x(), vy = particle.u.y(), vz = particle.u.z();
	double x = particle.p.x() - (x0_ + psx_ * particle.region.gidx_[0]);
	double y = particle.p.y() - (y0_ + psy_ * particle.region.gidx_[1]);
	double z = particle.p.z() - (z0_ + psz_ * particle.region.gidx_[2]);

	if (vx < 0) dx = -x / vx;
	else if (vx > 0) dx = (psx_ - x) / vx;

	if (vy < 0) dy = -y / vy;
	else if (vy > 0) dy = (psy_ - y) / vy;

	if (vz < 0) dz = -z / vz;
	else if (vz > 0) dz = (psz_ - z) / vz;

	if (dx < dy && dx < dz)
	{
		gidxNext[0] = particle.region.gidx_[0] + (vx < 0 ? -1 : 1);
		if (gidxNext[0] < 0 || gidxNext[0] >= (short)nx_) idx = -1;
		else idx = (gidxNext[2] * ny_ + gidxNext[1]) * nx_ + gidxNext[0];
		return dx;
	}
	else if (dy < dx && dy < dz)
	{
		gidxNext[1] = particle.region.gidx_[1] + (vy < 0 ? -1 : 1);
		if (gidxNext[1] < 0 || gidxNext[1] >= (short)ny_) idx = -1;
		else idx = (gidxNext[2] * ny_ + gidxNext[1]) * nx_ + gidxNext[0];
		return dy;
	}
	else
	{
		gidxNext[2] = particle.region.gidx_[2] + (vz < 0 ? -1 : 1);
		if (gidxNext[2] < 0 || gidxNext[2] >= (short)nz_) idx = -1;
		else idx = (gidxNext[2] * ny_ + gidxNext[1]) * nx_ + gidxNext[0];
		return dz;
	}
}

double mcPTGrid3D::getDNearInsideVoxel(const mcParticle& particle)
{
	double dx = particle.p.x() - x0_ - psx_ * particle.region.gidx_[0];
	double dy = particle.p.y() - y0_ - psy_ * particle.region.gidx_[1];
	double dz = particle.p.z() - z0_ - psz_ * particle.region.gidx_[2];
	if (dx < -GEOM_EPSILON || dx > psx_ + GEOM_EPSILON || dy < -GEOM_EPSILON || dy > psy_ + GEOM_EPSILON || dz < -GEOM_EPSILON || dz > psz_ + GEOM_EPSILON)
		throw std::exception("mcPTGrid3D::getDNearInsideVoxel: particle outside expected volume");
	return fabs(min(min(min(dx, psx_ - dx), min(dy, psy_ - dy)), min(dz, psz_ - dz)));
}

void mcPTGrid3D::beginTransport(mcParticle& p)
{
	// Указатель на транспортный объект, в котором частица находится в данный момент
	p.transport_ = this;

	mcParticle* particle = p.thread_->NextParticle();
	*particle = p;
	particle->p = p.p * mwtot_;
	particle->plast = p.plast * mwtot_;
	particle->u = particle->u.transformDirection(mwtot_);
	particle->dnear = 0;
	particle->mfps = HowManyMFPs(p.thread_->rng());

	// Переместить частицу на поверхность, если она еще не внутри
	int idx = getIdxAtPoint(particle->p, particle->region.gidx_);
	if (idx < 0)
	{
		double f = getDistanceOutside(*particle);
		if (f == DBL_MAX)
			return endTransport(particle);
		particle->p += particle->u * (f + FLT_EPSILON);
		idx = getIdxAtPoint(particle->p, particle->region.gidx_);
	}
	particle->region.idx_ = idx;
	particle->region.medidx_ = (short)mcbody_.GetMediaIdx(idx);
	particle->regDensityRatio = mcbody_.GetDensity(idx);

	// Транспорт в локальной системе координат
	simulate(p.thread_);
}

void mcPTGrid3D::beginTransportInside(mcParticle& p)
{
	beginTransport(p);
}

mc_move_result_t mcPTGrid3D::moveParticle(mcParticle* particle, double& step, double& edep)
{
	edep = 0;
	step = 0;

	// Важный момент! 
	// Если частица находится за пределами (на что указывает индекс региона), то возвращаемс с метокй о выходе.
	if (particle->region.idx_ < 0)
	{
		particle->exitSurface_ = mcParticle::temb_shit_t::External;
		return MCMR_EXIT;
	}

	// HACK!!
	// По непонятным причинам координаты частицы могут быть абсурдными.
	// Удалаяем такие частицы
	if (_isnan(particle->p.x()) != 0)
	{
		//cout << "Non number position or direction in object: " << this->getName() << endl;
		cout << "Non number position in object: " << this->getName() << endl;
		cout << "Position: " << particle->p;
		cout << "Direction: " << particle->u;
		particle->thread_->RemoveParticle();
		return MCMR_DISCARGE;
	}

	const mcPhysics* phys = media_->getPhysics(particle->t);

	//Параметры сред транспорта фотонов и электронов
	const mcMedium* med = media_->getMedium(particle->t, particle->region.medidx_);

	// Частицы с энергией ниже критической должны быть уничтожены раньше любых расчетов транспорта.
	if (phys->Discarge(particle, *med, edep))
		return MCMR_DISCARGE;

	// Hack!!! GG 20171030
	if (_isnan(particle->ke) != 0)
	{
		cout << "Non number energy: " << this->getName() << endl;
		cout << "Position: " << particle->p;
		cout << "Direction: " << particle->u;
		particle->thread_->RemoveParticle();
		return MCMR_DISCARGE;
	}

	double freepath = phys->MeanFreePath(particle->ke, *med, particle->regDensityRatio);
	step = freepath * particle->mfps;

	// Расстояние до границы ячейки в направлении частицы.
	// Если потом мы сделаем шаг до границы, то должны переместить частиу в следующую ячейку.
	// Именно в это момент мы пометим частицу как находящуюся на грани
	short gidxNext[3];
	int idx;
	double dist = getDistanceInsideVoxel(*particle, gidxNext, idx) + GEOM_EPSILON;

	// Для заряженных частиц dnear критично, так как код может убивать  частицу,
	// если ее range меньше расстояния до границы воксела, что верно в большинстве случаев.
	if (particle->t != mc_particle_t::MCP_PHOTON)
		particle->dnear = getDNearInsideVoxel(*particle);

	if (step < dist)
	{
		double stepRequested = step;
		edep = phys->TakeOneStep(particle, *med, step);

		if (step < stepRequested)
			return MCMR_CONTINUE;
		else
			return MCMR_INTERUCT;
	}
	else
	{
		// HACK! На поверхности возможно залипание, если расстояние в пределах погрешности вычислений.
		if (dist < GEOM_EPSILON)
			dist = GEOM_EPSILON;
		step = dist;
		edep = phys->TakeOneStep(particle, *med, step);

		particle->mfps -= step / freepath;
		if (step == dist)
		{
			// Добрались до границы воксела. Нужно перенастороить частицу на новый.
			particle->region.idx_ = idx;
			if (idx >= 0)
			{
				particle->region.medidx_ = (short)mcbody_.GetMediaIdx(idx);
				particle->regDensityRatio = mcbody_.GetDensity(idx);
				particle->region.gidx_[0] = gidxNext[0];
				particle->region.gidx_[1] = gidxNext[1];
				particle->region.gidx_[2] = gidxNext[2];
			}
		}
		return MCMR_CONTINUE;
	}
}

std::shared_ptr<std::vector<double>> mcPTGrid3D::GetSummDose() const
{
	auto M = std::make_shared<std::vector<double>>(size_, 0);
	// Сумма энергий
	for (auto& v : ((mcScorePTGrid3D*)score_)->EnergyMatrix())
	{
		for (unsigned i = 0; i < size_; i++)
			(*M)[i] += v[i];
	}

	// Деление на массу
	auto volume = psx_ * psy_ * psz_;	// cm^3
	for (unsigned i = 0; i < size_; i++)
	{
		double fd = mcbody_.GetPhysicalDensity(i);

		// Обнуляем дозы в воздухе
		if (fd < 0.01)
			(*M)[i] = 0;
		else
			(*M)[i] /= volume * fd;
	}

	return M;
}

void mcPTGrid3D::ClearEnergyMatrix()
{
	((mcScorePTGrid3D*)score_)->ClearEnergyMatrix();
}

void mcPTGrid3D::dump(ostream& os) const
{
	__super::dump(os);
}
