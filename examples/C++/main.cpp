
#include <SHARPlib/constants.h>
#include <SHARPlib/layer.h>
#include <SHARPlib/parcel.h>
#include <SHARPlib/thermo.h>
#include <SHARPlib/winds.h>

#include <chrono>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

// This is just a simple struct for passing arroun
// array pointers that exists for convenience. It is an
// example of how one could make their own data structure
// to retain state, but may not be the best approach for
// your particular use case...
//
// We use uniwue pointers because they
// are managed memory to raw array pointers,
// and the library takes raw array pointers.
// You could also use a std::vector and pass
// in the .data() call, which returns the
// raw pointer
struct Profile {
    std::unique_ptr<float[]> pres;
    std::unique_ptr<float[]> hght;
    std::unique_ptr<float[]> tmpk;
    std::unique_ptr<float[]> dwpk;
    std::unique_ptr<float[]> mixr;
    std::unique_ptr<float[]> vtmp;
    std::unique_ptr<float[]> theta;
    std::unique_ptr<float[]> theta_e;
    std::unique_ptr<float[]> uwin;
    std::unique_ptr<float[]> vwin;
    std::unique_ptr<float[]> wspd;
    std::unique_ptr<float[]> wdir;

    std::unique_ptr<float[]> pcl_vtmp;
    std::unique_ptr<float[]> pcl_buoy;
    const size_t NZ;
    Profile(size_t NZ) : NZ(NZ) {
        pres = std::make_unique<float[]>(NZ);
        hght = std::make_unique<float[]>(NZ);
        tmpk = std::make_unique<float[]>(NZ);
        dwpk = std::make_unique<float[]>(NZ);
        mixr = std::make_unique<float[]>(NZ);
        vtmp = std::make_unique<float[]>(NZ);
        theta = std::make_unique<float[]>(NZ);
        theta_e = std::make_unique<float[]>(NZ);
        uwin = std::make_unique<float[]>(NZ);
        vwin = std::make_unique<float[]>(NZ);
        wspd = std::make_unique<float[]>(NZ);
        wdir = std::make_unique<float[]>(NZ);

        // To compute CAPE/CINH, we need some scratch
        // arrays to store parcel virtual temperature
        // and buoyancy
        pcl_vtmp = std::make_unique<float[]>(NZ);
        pcl_buoy = std::make_unique<float[]>(NZ);
    };
};

// splits a string on a delimiter and stores in a vector
std::vector<std::string> split(std::string& s, std::string delimiter) {
    size_t pos_start = 0;
    size_t pos_end = 0;
    size_t delim_len = delimiter.length();
    std::string token;
    std::vector<std::string> res;

    while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
        token = s.substr(pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back(token);
    }

    res.push_back(s.substr(pos_start));
    return res;
}

void build_profile(Profile& prof, std::vector<std::string>& row, int idx) {
    float pres = std::stof(row[0]) * sharp::HPA_TO_PA;
    float hght = std::stof(row[1]);
    float tmpk = std::stof(row[2]) + sharp::ZEROCNK;
    float dwpk = std::stof(row[3]) + sharp::ZEROCNK;
    float wdir = std::stof(row[4]);
    float wspd = std::stof(row[5]);

    prof.pres[idx] = pres;
    prof.hght[idx] = hght;
    prof.tmpk[idx] = tmpk;
    prof.dwpk[idx] = dwpk;
    prof.wdir[idx] = wdir;
    prof.wspd[idx] = wspd;

    prof.mixr[idx] = sharp::mixratio(pres, dwpk);
    prof.vtmp[idx] = sharp::virtual_temperature(tmpk, prof.mixr[idx]);
    prof.theta[idx] = sharp::theta(pres, tmpk, sharp::THETA_REF_PRESSURE);
    prof.theta_e[idx] = sharp::thetae(pres, tmpk, dwpk);

    sharp::WindComponents uv = sharp::vector_to_components(wspd, wdir);

    prof.uwin[idx] = uv.u;
    prof.vwin[idx] = uv.v;
}

void read_sounding(std::unique_ptr<Profile>& prof, std::string filename) {
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

    /*Profile* prof = new Profile(NLINES);*/
    prof = std::make_unique<Profile>(NLINES);

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

                build_profile(*prof, row, idx);

                idx += 1;
            }
            // see if the data section begind
            if (line == begin) found_begin = true;
        }
        sndfile.close();
        std::cout << "Success reading: " << filename << std::endl;
        std::cout << "Number of vertical levels: " << prof->NZ << std::endl;

    }

    else {
        std::cout << "Unable to open file: " << filename << std::endl;
    }
}

int main() {
    std::string snd_file1 =
        "../../data/test_snds/20160524_2302_EF3_37.57_-100.13_108_613967.snd";
    std::string snd_file2 = "../../data/test_snds/hires-SPC.txt";
    std::unique_ptr<Profile> prof;
    read_sounding(prof, snd_file1);

    if (prof) {
        static sharp::lifter_cm1 lifter;
        lifter.ma_type = sharp::adiabat::adiab_ice;
        auto start_time = std::chrono::system_clock::now();

        // This is the PressureLayer over which to search
        // for the Most Unstable parcel
        sharp::PressureLayer mu_lyr = {prof->pres[0], prof->pres[0] - 40000.0f};

        // This is the PressureLayer used to define the mixing
        // depth for the mixed-layer parcel
        sharp::PressureLayer ml_lyr = {prof->pres[0], prof->pres[0] - 10000.0f};

        sharp::Parcel sfc_pcl = sharp::Parcel::surface_parcel(
            prof->pres[0], prof->tmpk[0], prof->dwpk[0]);

        // The MU parcel already has its CAPE/CINH set since
        // it has to search based on CAPE anyway.
        sharp::Parcel mu_pcl = sharp::Parcel::most_unstable_parcel(
            mu_lyr, lifter, prof->pres.get(), prof->hght.get(),
            prof->tmpk.get(), prof->vtmp.get(), prof->dwpk.get(),
            prof->pcl_vtmp.get(), prof->pcl_buoy.get(), prof->NZ);

        sharp::Parcel ml_pcl = sharp::Parcel::mixed_layer_parcel(
            ml_lyr, prof->pres.get(), prof->hght.get(), prof->theta.get(),
            prof->mixr.get(), prof->NZ);

        // lift the parcel to get its virtual temperature,
        // compute the buoyancy, and then integrate to get
        // the CAPE/CINH
        sfc_pcl.lift_parcel(lifter, prof->pres.get(), prof->pcl_vtmp.get(),
                            prof->NZ);
        sharp::buoyancy(prof->pcl_vtmp.get(), prof->vtmp.get(),
                        prof->pcl_buoy.get(), prof->NZ);
        sfc_pcl.cape_cinh(prof->pres.get(), prof->hght.get(),
                          prof->pcl_buoy.get(), prof->NZ);

        ml_pcl.lift_parcel(lifter, prof->pres.get(), prof->pcl_vtmp.get(),
                           prof->NZ);
        sharp::buoyancy(prof->pcl_vtmp.get(), prof->vtmp.get(),
                        prof->pcl_buoy.get(), prof->NZ);
        ml_pcl.cape_cinh(prof->pres.get(), prof->hght.get(),
                         prof->pcl_buoy.get(), prof->NZ);

        auto end_time = std::chrono::system_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
                            end_time - start_time)
                            .count();

        std::cout << "Lifting 3 parcels took: " << duration << "us"
                  << std::endl;

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
}
