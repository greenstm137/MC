// Radiation Oncology Monte Carlo open source project
//
// Author: [2019] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include <vector>
#include <memory>

enum media_type_e {
	MT_Air = 0,
	MT_Water = 1,
	MT_Brain = 2,
	MT_Bone = 3,
	MT_Skin = 4
};

struct IsotopeElement
{
	// јтомный вес изотопа
	double A;

	// јтомный номер (зар€д) изотопа
	unsigned short Z;

	//  оличество атомов в среде на один элемент материала.
	// ¬ описании сред соблюдаютс€ относительные пропорции
	// атомов, но нормировки может не быть.
	// Ёто следует учитывать при расчете сечений взаимодействий дл€ сред.
	double W;
};

struct MediaComposition
{
	// »дентификатор замен€ет текстовое название, с которым проблема передачи текстов 
	// через границы в том числе и из-за разных представлений текстов
	media_type_e Id;

	std::vector<IsotopeElement> Elements;

	// —тандартна€ плотность материала среды (g/cm^3)
	double Dens = 0;
};

//  ласс представлени€ пиксельного фантома (человека) .
class mcPTBody
{
public:
	mcPTBody();
	virtual ~mcPTBody(void);

	void SetGrid(unsigned nx, unsigned ny, unsigned nz,
		double x0, double y0, double z0, double psx, double psy, double psz);

	// ћассив индексов сред
	std::vector<unsigned short>& GetMediaIndexes() { return *mediaIdxs_; }

	// ћассив плотностей относительно плотности среды
	std::vector<double>& GetDensities() { return *dens_; }

	unsigned short GetMediaIdx(int i) const { return mediaIdxs_->at(i); }
	double GetDensity(int i) const { return dens_->at(i); }
	double GetPhysicalDensity(int i) const { return dens_->at(i) * defaultDensities_->at(mediaIdxs_->at(i)); }

	// ”становка стандартных плотностей материалов из PEGS4
	void SetDefaultDencities(const std::vector<double>& densities);

public:
	unsigned nx_, ny_, nz_;
	unsigned size_;
	double psx_, psy_, psz_;
	// first voxel coordinates (others should have greater or equal)
	double x0_, y0_, z0_;

	std::unique_ptr<std::vector<unsigned short>> mediaIdxs_;
	std::unique_ptr<std::vector<double>> dens_;
	std::unique_ptr<std::vector<double>> defaultDensities_;
};
