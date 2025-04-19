// Radiation Oncology Monte Carlo open source project
//
// Author: [2020] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcscore.h"
#include <vector>

// Класс скоринга для поддержки транспорта в 3D сетке типа анатомии человека.
// Скоринг вынесен из объекта транспорта, хотя и напрашивался туда,
// для удовлетворения стандартному интерфейсу.
class mcScorePTGrid3D : public mcScore
{
public:
	mcScorePTGrid3D(const char* module_name, int nThreads, int nx, int ny, int nz);
	virtual ~mcScorePTGrid3D(void);

	void ScoreFluence(const mcParticle& particle) override;

	void ScorePoint(double edep
		, int iThread
		, const mcRegionReference& region
		, mc_particle_t pt
		, const geomVector3D& p0) override;

	void ScoreLine(double edep
		, int iThread
		, const mcRegionReference& region
		, mc_particle_t pt
		, const geomVector3D& p0
		, const geomVector3D& p1) override;

	void CE2D() override;

	void dumpStatistic(ostream&) const override;

	void ClearEnergyMatrix();
	const std::vector<std::vector<double>>& EnergyMatrix() const;

protected:
	unsigned nx_, ny_, nz_;
	// Массив накопленных энергий
	std::unique_ptr<std::vector<std::vector<double>>> energy_;
};
