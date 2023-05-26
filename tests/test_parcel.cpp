#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>

#include "doctest.h"
#include <SHARPlib/constants.h>
#include <SHARPlib/profile.h>
#include <SHARPlib/winds.h>
#include <SHARPlib/parcel.h>

// splits a string on a delimiter and stores in a vector
std::vector<std::string> split(std::string& s, std::string delimiter) {
    size_t pos_start = 0;
    size_t pos_end = 0;
    size_t delim_len = delimiter.length();
    std::string token;
    std::vector<std::string> res;

    while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }

    res.push_back (s.substr (pos_start));
    return res;
}

void build_profile(sharp::Profile* prof, std::vector<std::string>& row, int idx) {
    float pres = std::stof(row[0]) * 100.0;
    float hght = std::stof(row[1]);
    float tmpk = std::stof(row[2]) + sharp::ZEROCNK;
    float dwpk = std::stof(row[3]) + sharp::ZEROCNK;
    float wdir = std::stof(row[4]);
    float wspd = std::stof(row[5]);

    prof->pres[idx] = pres; 
    prof->hght[idx] = hght;
    prof->tmpk[idx] = tmpk;
    prof->dwpk[idx] = dwpk;
    prof->wdir[idx] = wdir;
    prof->wspd[idx] = wspd;

    prof->vtmp[idx] = sharp::virtual_temperature(pres, tmpk, dwpk);
    prof->mixr[idx] = sharp::mixratio(pres, dwpk);
    prof->theta[idx] = sharp::theta(pres, tmpk, sharp::THETA_REF_PRESSURE);
    prof->theta_e[idx] = sharp::thetae(pres, tmpk, dwpk);
    
    sharp::WindComponents uv = sharp::vector_to_components(wspd, wdir);

    prof->uwin[idx] = uv.u;
    prof->vwin[idx] = uv.v;
}

sharp::Profile* read_sounding(std::string filename) {
    std::ifstream sndfile(filename);
    std::string line;

    std::string begin = "%RAW%";
    std::string end = "%END%";
    bool found_begin = false;
    bool found_end = false;
    int NLINES = 0;
    // first iteration - count the number of data
    // rows so that we can allocate arrays of appropraite
    // size to store the data
    if (sndfile.is_open()) {
        while (std::getline(sndfile, line)) {
            if (line == end) found_end = true;

            // check if we are within the data section
            if (found_begin && !found_end) {
                NLINES += 1;
            }
            // see if the data section begind
            if (line == begin) found_begin = true;
        }
    }

    sharp::Profile* prof = new sharp::Profile(NLINES, sharp::Source::PFC);

    // return to the beginning of the file
    sndfile.clear();
    sndfile.seekg(0);

    found_begin = false;
    found_end = false;
    int idx = 0;
    // now loop again, this time 
    // splitting the data
    if (sndfile.is_open()) {
        // iterate over the lines in the file
        while (std::getline(sndfile, line)) {
            // loop-logic: if we're at the end of the
            // data section, set that first. 
            if (line == end) found_end = true;

            // check if we are within the data section
            if (found_begin && !found_end) {

                // split the line on the comma
                std::vector row = split(line, ",");

                build_profile(prof, row, idx);

                idx += 1;
            }
            // see if the data section begind
            if (line == begin) found_begin = true;
        }
        sndfile.close();
        std::cout << "Success reading: " << filename << std::endl;
        std::cout << "Number of vertical levels: " << prof->NZ << std::endl;
    
        return prof;
    }

    else {
        std::cout << "Unable to open file: " << filename << std::endl;
        return nullptr;
    }
}


TEST_CASE("Testing new parcel definitions") {
    std::string fname1 = "/users/khalbert/CODEBASE/NSHARP-server/unprocessed/sars_supercell/99050323f0.okc";
    std::string fname2 = "/users/khalbert/Downloads/newSPC.txt";
    std::string fname3 = "/users/khalbert/22112200.KEY";
	std::string fname4 = "/home/kthalbert/CODEBASE/SHARP-calc/swig/test/hires-SPC.txt";
	std::string fname5 = "/users/khalbert/CODEBASE/NSHARP-calc/swig/test/hires-SPC.txt";


    auto start_time = std::chrono::system_clock::now();
    sharp::Profile* prof = read_sounding(fname5);
    auto end_time = std::chrono::system_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();

    std::cout << "Reading file took: " << duration << "ms" << std::endl;
    if (prof) {
        sharp::Parcel sfc_pcl;
        sharp::Parcel mu_pcl;
        sharp::Parcel ml_pcl;

        sharp::define_parcel(prof, &sfc_pcl, sharp::LPL::SFC);
        sharp::define_parcel(prof, &mu_pcl, sharp::LPL::MU);
        sharp::define_parcel(prof, &ml_pcl, sharp::LPL::ML);

        start_time = std::chrono::system_clock::now();

        sharp::parcel_wobf(prof, &sfc_pcl);
        sharp::parcel_wobf(prof, &ml_pcl);
        sharp::parcel_wobf(prof, &mu_pcl);

        end_time = std::chrono::system_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();

        std::cout << "Lifting 3 parcels took: " << duration << "us" << std::endl;
        std::cout << "Average was: " << duration / 3.0 << "us" << std::endl;

        std::cout << "SFC PCL\t";
        std::cout << "LFC PRES: " << sfc_pcl.lfc_pressure << "\t"; 
        std::cout << "EL PRES: " << sfc_pcl.eql_pressure << "\t"; 
        std::cout << "CAPE: " << sfc_pcl.cape << "\t"; 
        std::cout << "CINH: " << sfc_pcl.cinh << std::endl;

        std::cout << "MU PCL\t";
        std::cout << "LFC PRES: " << mu_pcl.lfc_pressure << "\t"; 
        std::cout << "EL PRES: " << mu_pcl.eql_pressure << "\t"; 
        std::cout << "CAPE: " << mu_pcl.cape << "\t"; 
        std::cout << "CINH: " << mu_pcl.cinh << std::endl;

        std::cout << "ML PCL\t";
        std::cout << "LFC PRES: " << ml_pcl.lfc_pressure << "\t"; 
        std::cout << "EL PRES: " << ml_pcl.eql_pressure << "\t"; 
        std::cout << "CAPE: " << ml_pcl.cape << "\t"; 
        std::cout << "CINH: " << ml_pcl.cinh << std::endl;
    }
    delete prof;
}
