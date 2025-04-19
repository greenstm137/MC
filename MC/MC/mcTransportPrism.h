// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include "mctransport.h"

//  ласс транспорта в параллелепипеде.
// ÷ентр системы координат объекта находитс€
// в середине грани основани€ (Z=const)
class mcTransportPrism : public mcTransport
{
public:
	mcTransportPrism(void);
	mcTransportPrism(const geomVector3D& orgn, const geomVector3D& vz, const geomVector3D& vx,
		double ax, double ay, double az);
	mcTransportPrism(const geomVector3D& orgn, const geomVector3D& vz, const geomVector3D& vx);
	virtual ~mcTransportPrism(void);

	void setGeometry(double ax, double ay, double az);

	// ”становка в новой геометрии, более подход€щей дл€ параллелепипеда, содержащего пациента.
	void setGeometry(double ax, double ay, double az, double x0, double y0, double z0);

	double ax() const { return ax_; }
	double ay() const { return ay_; }
	double az() const { return az_; }

	double getDistanceInside(mcParticle& p) const override;
	double getDistanceOutside(mcParticle& p) const override;

	void dump(ostream& os) const override;
	void dumpVRML(ostream& os)const override;

protected:
	double getDNearInside(const geomVector3D& p) const override;

protected:
	double ax_;
	double ay_;
	double az_;

	//  оординаты нулевого угла в собственной системе относительно точки, на которую указывает orgn
	double x0_, y0_, z0_;

	// —мещение между старым и новым взгл€дом дл€ ускорени€ расчетов
	double dx0_, dy0_, dz0_;
};
