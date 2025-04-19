// Radiation Oncology Monte Carlo open source project
//
// Author: [2019-20] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mcTransportPrism.h"

class mcPTBody;
class mcScorePTGrid3D;

// Класс транспорта в воксельной 3D геометрии, представляющей в том числе и человека.
class mcPTGrid3D : public mcTransportPrism
{
public:
	mcPTGrid3D(const geomVector3D& orgn, const geomVector3D& vz, const geomVector3D& vx, const mcPTBody& mcbody, int nThreads);
	virtual ~mcPTGrid3D(void);

	// В отличие от стандартного начала переносит частицу на поверхность фантома 
	// и вычисляет индекс стартовой ячейки
	void beginTransport(mcParticle& p) override;

	// Старт транспорта может вызываться как в объекте, вложенном во внутреннюю структуру - воздушеый слой.
	void beginTransportInside(mcParticle& p) override;

	// Виртуальная функция перемещения частицы полностью покрывает специфику транспорта в сетке.
	// Предыдущее решение через собственную функцию транспорта бло полностью ошибочным
	// так как любой симулятор опусташал стэк всего потока и для частиц вперемешку вызывались разные 
	// симуляторы, в том числе и не соответствующие данному объекта транспорта.
	mc_move_result_t moveParticle(mcParticle* particle, double& step, double& edep) override;

	// Возвращает дозу, как сумму энергий по потокам, деленную на массу воксела
	std::shared_ptr<std::vector<double>> GetSummDose() const;

	void ClearEnergyMatrix();

	unsigned nx() const { return nx_; }
	unsigned ny() const { return ny_; }
	unsigned nz() const { return nz_; }

	double x0() const { return x0_; }
	double y0() const { return x0_; }
	double z0() const { return z0_; }

	double psx() const { return psx_; }
	double psy() const { return psy_; }
	double psz() const { return psz_; }

	void dump(ostream& os) const override;

protected:
	int getIdxAtPoint(const geomVector3D& p, short* pgidx) const;
	double getDistanceInsideVoxel(const mcParticle& particle, short* gidxNext, int& idx);
	double getDNearInsideVoxel(const mcParticle& particle);

	const mcPTBody& mcbody_;
	int nThreads_;

	unsigned nx_, ny_, nz_;
	unsigned size_;
	// voxel sizes
	double psx_, psy_, psz_;
	// first voxel coordinates
	double x0_, y0_, z0_;
};
