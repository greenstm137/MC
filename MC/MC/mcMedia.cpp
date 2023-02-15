#include "mcMedia.h"
#include "mcMediumXE.h"
#include "mcMediumProton.h"
#include "mcMediumNeutron.h"
#include "mcPhysicsPhoton.h"
#include "mcPhysicsElectron.h"
#include "mcPhysicsPositron.h"
#include "mcPhysicsProton.h"
#include "mcPhysicsNeutron.h"
#include "mcParticle.h"
#include "mcCSNuclear.h"
#include "../geometry/text.h"
#include <fstream>
//#include <string>
//#include <iostream>
#include <filesystem>

namespace fs = std::filesystem;

mcMedia::mcMedia(void)
{
	physics_.resize(MCP_NTYPES, nullptr);
	physics_[MCP_PHOTON] = new mcPhysicsPhoton();
	physics_[MCP_NEGATRON] = new mcPhysicsElectron();
	physics_[MCP_POSITRON] = new mcPhysicsPositron();
	physics_[MCP_PROTON] = new mcPhysicsProton();
	physics_[MCP_NEUTRON] = new mcPhysicsNeutron();
}

mcMedia::~mcMedia(void)
{
	int i;
	for (i = 0; i < (int)xes_.size(); i++)
		delete xes_[i];
	for (i = 0; i < (int)protons_.size(); i++)
		delete protons_[i];
	for (i = 0; i < (int)neutrons_.size(); i++)
		delete neutrons_[i];
	for (i = 0; i < (int)physics_.size(); i++)
		delete physics_[i];
}

void mcMedia::addName(const char* mname)
{
	if (!xes_.empty() || !protons_.empty() || !neutrons_.empty())
		throw std::exception("Can't add media names after data initialization");
	mnames_.push_back(mname);
}

short mcMedia::getMediumIdx(const char* mname) const
{
	for (short i = 0; i < (short)mnames_.size(); i++)
		if (mnames_[i] == mname) return i;
	throw std::exception((string("Medium \"") + mname + string("\" does not exist")).c_str());
}

const mcMediumXE* mcMedia::getMediumXE(short idx) const
{
	if (idx >= (short)xes_.size())
		throw std::exception("Medium index for photons an electrons is too big");
	return (mcMediumXE*)xes_[idx];
}

const mcMediumProton* mcMedia::getProtonMedium(short idx) const
{
	if (idx >= (short)protons_.size())
		throw std::exception("Medium index for protons is too big");
	return (mcMediumProton*)protons_[idx];
}

const mcMediumProton* mcMedia::getNeutronMedium(short idx) const
{
	if (idx >= (short)neutrons_.size())
		throw std::exception("Medium index for protons is too big");
	return (mcMediumProton*)neutrons_[idx];
}

void mcMedia::initXEFromStream(istream& is)
{
	if (!xes_.empty())
		throw std::exception("Photon and electron crossections already initialized");
	for (int i = 0; i < (int)mnames_.size(); i++)
		xes_.push_back(new mcMediumXE());

	// ������ ������
	string line, s1, s2;
	getline(is, line, '\n');
	while (!is.fail())
	{
		if (line.find("MEDIUM=") != string::npos)
		{
			GetTwoStringsFromLine(line, s1, s2);
			GetTwoStringsFromLine(s2, line, s1);

			// ���������, ����� �� ������ ����� ��� ��������?
			int i;
			for (i = 0; i < (int)mnames_.size(); i++)
				if (mnames_[i] == line) break;

			if (i < (int)mnames_.size()) {
				xes_[i]->name_ = line;
				((mcMediumXE*)xes_[i])->read(is);
			}
		}
		getline(is, line, '\n');
	}

	// ���������, ��� �� ����� ���������
	string errmedia;
	for (int i = 0; i < (int)xes_.size(); i++)
	{
		if (xes_[i]->status_ != mcMedium::LOADED) {
			errmedia += xes_[i]->name_;
			errmedia += "\n";
		}
	}
	if (!errmedia.empty())
		throw std::exception((string("The following XE media were not loaded succcessfuly:\n") + errmedia).c_str());
}

void mcMedia::initXEFromFile(const string& fname)
{
	ifstream is(fname.c_str());
	if (is.fail())
		throw std::exception((string("Can't open XE data file: ") + fname).c_str());
	initXEFromStream(is);
}

void mcMedia::initProtonDeDxFromStream(istream& is)
{
	if (!protons_.empty())
		throw std::exception("Proton crossectons already initialized");
	int i;
	for (i = 0; i < (int)mnames_.size(); i++)
		protons_.push_back(new mcMediumProton());

	// ������ ������ - ����� � ���� ������� ��������� ���������� XA, ������ ��������� �������� ������
	string line, s1, s2, s3, s4;
	getline(is, line, '\n');
	while (!is.fail())
	{
		if (line.find("MEDIUM=") != string::npos)
		{
			GetTwoStringsFromLine(line, s1, s2);
			GetTwoStringsFromLine(s2, line, s1);

			// ���������, ����� �� ������ ����� ��� ��������?
			int i;
			for (i = 0; i < (int)mnames_.size(); i++)
				if (mnames_[i] == line) break;

			if (i < (int)mnames_.size()) {
				// ������������� ��������� ������ input file VER=0.0.0
				GetTwoStringsFromLine(s1, s2, s3);
				GetTwoStringsFromLine(s3, s1, s4);
				if ((s2 == "VER") || (s3 == "0.0.0")) {
					protons_[i]->name_ = line;
					((mcMediumProton*)protons_[i])->read(is);
				}
				else {
					//throw std::exception("Wrong Proton media data version"); 
					//� �������� ������ ����� ���� ������ � ���� �� ����� � ������ ������, 
					// ��� ��� ������ �� ��������� ������
				}
			}
		}
		getline(is, line, '\n');
	}

	// ���������, ��� �� ����� ���������
	string errmedia;
	for (int i = 0; i < (int)protons_.size(); i++)
	{
		if (protons_[i]->status_ != mcMedium::LOADED) {
			errmedia += mnames_[i];
			errmedia += "\n";
		}
	}
	if (!errmedia.empty())
		throw std::exception((string("The following Proton media were not loaded succcessfuly:\n") + errmedia).c_str());
}

void mcMedia::initProtonFromFiles(const string& fname, const string& pstardir, const string& icru63dir)
{
	// ������ ������� ��������� ������������� (����������, 2008)
	ifstream is(fname.c_str());
	if (is.fail())
		throw std::exception((string("Can't open Proton data file: ") + fname).c_str());
	initProtonDeDxFromStream(is);

	// TODO: ����� ������ ��������� ������������ �� ���� ������ PSTAR.
	// Note: ��� ���� �������� ������ ������������ ����� ������, 
	// �� ������� ���� � ������ � ��������� �� ����� ����� � ���������.
	// ����� ����������, ��� �� ���-�� ������ ��� ���� ������ ��� ��� ������ �������
	// ��� ������� ����� ������ ����� ���������� ��������� I.
	// 
	// �� �����������! ������ �������� ������ �������.
	//

	// ������� ��������� ���� ������ ������� ���� ��������� ������.
	// ����� ��� �������� ������������� ����� "P*.DAT".

	// ���� �� ���� ������ � ��������� ���� ������ ICRU 63 ��� ��������.

	//fs::path p = fs::current_path();
	//std::cout << "The current path " << p << " decomposes into:\n"
	//	<< "root-path " << p.root_path() << '\n'
	//	<< "relative path " << p.relative_path() << '\n';

	//std::filesystem::current_path("C:/MCSimulations/protons/");

	/*
	std::vector<FileEntry> CollectFiles(const fs::path & inPath)
	{
		std::vector<fs::path> paths;
		if (fs::exists(inPath) && fs::is_directory(inPath))
		{
			std::filesystem::recursive_directory_iterator dirpos{ inPath };

			std::copy_if(begin(dirpos), end(dirpos), std::back_inserter(paths),
				[](const fs::directory_entry& entry) {
					return entry.is_regular_file();
				}
			);
		}
		std::vector<FileEntry> files(paths.size());
		std::transform(paths.cbegin(), paths.cend(), files.begin(), FileEntry::Create);
		return files;
	}
	*/

	// Separators
	std::string energySeparator("+++++++++++++++++++++++++++++++");
	std::string evalueSeparator("Cross section summary information for");
	std::string tcsvalueSeparator("Nonelastic cross section =");
	std::string protonsSeparator("*protons*");
	std::string neutronsSeparator("*neutrons*");
	std::string protonCrossSectionSeparator("proton production cross section (mb) =");
	std::string neutronCrossSectionSeparator("neutron production cross section (mb) =");

	// ������, � ������� ������� ���������� ��� ��� ������ �������
	//auto dbData = std::make_unique<std::vector<std::unique_ptr<mcCSNuclear>>>();
	auto dbData = std::make_unique<std::vector<mcCSNuclear>>();

	// ���� �� ������ �������, � ������ �� ������� ���������� ������ ������ ��� ������ �������
	for (const auto& entry : fs::directory_iterator(icru63dir))
	{
		if (!fs::path(entry.path()).has_stem() || !fs::path(entry.path()).has_extension())
			continue;

		//string fname = fs::path(entry.path()).filename().string();
		string fname = fs::path(entry.path()).stem().string();
		string ext = fs::path(entry.path()).extension().string();

		if(fname[0] != 'P' || ext != ".DAT")
			continue;

		//if (std::regex_match(entry., reg))
		//std::cout << entry.path() << std::endl;
		//std::cout << elementName << std::endl;

		// ������ ���� ������� ��������.
		// ����� �������� �������� ����� �� ����� �����.
		string elementName = std::string(&fname[1]);

		ifstream isIcru(entry.path().c_str());
		if (isIcru.fail())
			throw std::exception((string("Can't open Proton data file: ") + entry.path().string()).c_str());

		// ���� ������ �������
		mcCSNuclear csForElement;
		csForElement.ElementName = elementName;

		// ������� ��� ������� �������� �������
		mcCSNuclearParticleEnergy csForEnergy;
		mcCSNuclearForAngleSpectrum particleAngle;
		mcCSNuclearForAngleSpectrum neutronAngle;

		// ������ ������ ������ ���� �� ������ � �������� ������ ����������
		string line, s1, s2, s3, s4;
		std::getline(isIcru, line, '\n');

		// ��������� �������� � ����� ����� �������� �� ��������� � ������ ��� ������������� ������
		bool isNewEnergyPrepare = false;
		bool isProtons = false;
		bool isNeutrons = false;
		int angleCount = 0;

		while (!isIcru.fail())
		{
			// ������ ����� �������
			if (line.find(energySeparator) != string::npos)
			{
				csForEnergy.Clear();
				isNewEnergyPrepare = true;
			}
			else if (isNewEnergyPrepare)
			{
				if (line.find(evalueSeparator) != string::npos)
				{
					csForEnergy.Energy = atof(&line[evalueSeparator.size()]);
				}
				else if (line.find(tcsvalueSeparator) != string::npos)
				{
					csForEnergy.TotalCrossSection = atof(&line[tcsvalueSeparator.size()]);
					isNewEnergyPrepare = false;
				}
			}

			// proton spectra
			else if (line.find(protonsSeparator) != string::npos)
			{
				csForEnergy.ProtonAngles.clear();
				isProtons = true;
				angleCount = 0;
			}
			else if (isProtons)
			{
				if (line.find("ANGLE(deg)") != string::npos)
				{
					std::vector<std::string> ss;
					GetStringArray(line, ss, " ");
					angleCount = ss.size() - 2;
					csForEnergy.ProtonAngles.resize(angleCount);
					for (int i = 0; i < angleCount; i++)
						csForEnergy.ProtonAngles[i].Angle = atof(ss[i + 1].c_str());
				}
				else if (line.find("ENERGY") != string::npos || line.find("------") != string::npos)
				{
				}
				else if (line.find(protonCrossSectionSeparator) != string::npos)
				{
					std::vector<std::string> ss;
					GetStringArray(line, ss, "=");
					csForEnergy.ProtonCrossSection = atof(ss[1].c_str());
					isProtons = false;
				}
				else
				{
					std::vector<std::string> ss;
					GetStringArray(line, ss, " ");
					for (int i = 0; i < (int)csForEnergy.ProtonAngles.size(); i++)
					{
						csForEnergy.ProtonAngles[i].SourceBinUps.push_back(atof(ss[1].c_str()));
						csForEnergy.ProtonAngles[i].SourceSpectrum.push_back(atof(ss[i + 2].c_str()));
					}
				}

			}

			// neutron spectra
			else if (line.find(neutronsSeparator) != string::npos)
			{
				csForEnergy.NeutronAngles.clear();
				isNeutrons = true;
				angleCount = 0;
			}
			else if (isNeutrons)
			{
				if (line.find("ANGLE(deg)") != string::npos)
				{
					std::vector<std::string> ss;
					GetStringArray(line, ss, " ");
					angleCount = ss.size() - 2;
					csForEnergy.NeutronAngles.resize(angleCount);
					for (int i = 0; i < angleCount; i++)
						csForEnergy.NeutronAngles[i].Angle = atof(ss[i + 1].c_str());
				}
				else if (line.find("ENERGY") != string::npos || line.find("------") != string::npos)
				{
				}
				else if (line.find(neutronCrossSectionSeparator) != string::npos)
				{
					std::vector<std::string> ss;
					GetStringArray(line, ss, "=");
					csForEnergy.NeutronCrossSection = atof(ss[1].c_str());
					isNeutrons = false;
				}
				else
				{
					std::vector<std::string> ss;
					GetStringArray(line, ss, " ");
					for (int i = 0; i < (int)csForEnergy.NeutronAngles.size(); i++)
					{
						csForEnergy.NeutronAngles[i].SourceBinUps.push_back(atof(ss[1].c_str()));
						csForEnergy.NeutronAngles[i].SourceSpectrum.push_back(atof(ss[i + 2].c_str()));
					}
				}
			}

			else if (line.find("**************************************************************************************************************************") != string::npos)
			{
				csForElement.Energies.push_back(csForEnergy);
			}
			std::getline(isIcru, line, '\n');
		}

		dbData->push_back(csForElement);
	}
}

void mcMedia::initNeutronFromStream(istream& is)
{
	if (!neutrons_.empty())
		throw std::exception("Neutron crossectons already initialized");
	int i;
	for (i = 0; i < (int)mnames_.size(); i++)
		neutrons_.push_back(new mcMediumNeutron());

	// ������ ������ - ����� � ���� ������� ��������� ���������� XA, ������ ��������� �������� ������
	string line, s1, s2, s3, s4;
	getline(is, line, '\n');
	while (!is.fail())
	{
		if (line.find("MEDIUM=") != string::npos)
		{
			GetTwoStringsFromLine(line, s1, s2);
			GetTwoStringsFromLine(s2, line, s1);

			// ���������, ����� �� ������ ����� ��� ��������?
			int i;
			for (i = 0; i < (int)mnames_.size(); i++)
				if (mnames_[i] == line) break;

			if (i < (int)mnames_.size()) {
				// ������������� ��������� ������ input file VER=0.0.0
				GetTwoStringsFromLine(s1, s2, s3);
				GetTwoStringsFromLine(s3, s1, s4);
				if ((s2 == "VER") || (s3 == "0.0.0")) {
					neutrons_[i]->name_ = line;
					((mcMediumNeutron*)neutrons_[i])->read(is);
				}
				else {
					//throw std::exception("Wrong Neutron media data version"); 
					//� �������� ������ ����� ���� ������ � ���� �� ����� � ������ ������, 
					// ��� ��� ������ �� ��������� ������
				}
			}
		}
		getline(is, line, '\n');
	}

	// ���������, ��� �� ����� ���������
	string errmedia;
	for (int i = 0; i < (int)neutrons_.size(); i++)
	{
		if (neutrons_[i]->status_ != mcMedium::LOADED) {
			errmedia += mnames_[i];
			errmedia += "\n";
		}
	}
	if (!errmedia.empty())
		throw std::exception((string("The following Neutron media were not loaded succcessfuly:\n") + errmedia).c_str());
}

void mcMedia::initNeutronFromFile(const string& fname)
{
	ifstream is(fname.c_str());
	if (is.fail())
		throw std::exception((string("Can't open Neutron data file: ") + fname).c_str());
	initNeutronFromStream(is);
}

const mcPhysics* mcMedia::getPhysics(int ptype) const
{
	if (ptype >= (int)physics_.size())
		throw std::exception("Unsupported particle type");
	return physics_[ptype];
}

const mcMedium* mcMedia::getMedium(int ptype, int idx) const
{
	if (idx >= (int)mnames_.size())
		throw std::exception("Unsupported particle type");
	if (ptype == MCP_PROTON)
		return protons_[idx];
	else if (ptype == MCP_NEUTRON)
		return neutrons_[idx];
	else
		return xes_[idx];
}
