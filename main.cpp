#include <iostream>
#include <CLI/CLI.hpp>
#include <spdlog/spdlog.h>
#include <hts.h>
#include <vcf.h>
#include <vector>


int main(int argc, char **argv)
{

    CLI::App app("Frequencies reader");
    std::string vcf_file_name;
    std::string out_file;
    app.add_option("-i,--input", vcf_file_name, "Input VCF file")->required();
    app.add_option("-o,--output", out_file, "Output file")->required();
    CLI11_PARSE(app, argc, argv);

    htsFile* inf = bcf_open(vcf_file_name.c_str(), "r");
    if (inf == NULL)
    {
        spdlog::error("Can't open vcf file: {}", vcf_file_name);
        std::exit(EXIT_FAILURE);
    }

    spdlog::info("Parsing vcf: {}", vcf_file_name);
 
    bcf_hdr_t* hdr = bcf_hdr_read(inf);

    bcf1_t* rec = bcf_init();
    if (rec == NULL)
    {
        spdlog::error("Error while parsing vcf file: {}", vcf_file_name);
        bcf_close(inf);
        bcf_hdr_destroy(hdr);
        std::exit(EXIT_FAILURE);
    }


    std::ofstream out_file_stream(out_file);
    if (not out_file_stream.is_open()) { spdlog::error("Can't open output file {}", out_file); }
    else { spdlog::info("Writing output to {}", out_file); }
    
    while (bcf_read(inf, hdr, rec) == 0)
    {
	std::string line;
	bcf_unpack(rec, BCF_UN_ALL);
    	int chrom = rec->rid;
        int position = rec->pos;
        const char* chr_name = bcf_hdr_id2name(hdr, rec->rid);
	bcf_dec_t qhat = rec ->d;
        char* ID = qhat.id;
        std::string someString(ID);
        out_file_stream << chr_name << "\t"<< position << "\t"<< ID << "\t";
	//get reference and alternate allele 
	std::string ref_allele = rec->d.allele[0];
	std::string alt_allele = rec->d.allele[1];
        bcf_info_t* info_pointer  = bcf_get_info(hdr, rec, "vep");
        std::string vep_val((char*)info_pointer->vptr);
	out_file_stream << ref_allele << "\t" << alt_allele << "\t";

	bcf_unpack(rec, BCF_UN_ALL);

	//AMR_MALE-------------------------------
	bcf_info_t* AN_amr_male  = bcf_get_info(hdr, rec, "AN_amr_male");
        bcf_info_t* AC_amr_male  = bcf_get_info(hdr, rec, "AC_amr_male");
        bcf_info_t* AF_amr_male  = bcf_get_info(hdr, rec, "AF_amr_male");
        int AC_amr_male2 = AC_amr_male->v1.i;
        int AN_amr_male2 = AN_amr_male->v1.i;
        float AF_amr_male2;
        if(AF_amr_male == NULL) {
            if(AN_amr_male2 != 0){
                AF_amr_male2 = AC_amr_male2 / AN_amr_male2;
            }
            else{
                AF_amr_male2 = 0;
            }
        }else{
            AF_amr_male2 = AF_amr_male->v1.f;
        }

        out_file_stream << "AC_amr_male=" << AC_amr_male2 <<";";
        out_file_stream << "AN_amr_male=" << AN_amr_male2 <<";";
        out_file_stream << "AF_amr_male=" << AF_amr_male2 << ";";

	
	//AMR_FEMALE
	 bcf_info_t* AN_amr_female  = bcf_get_info(hdr, rec, "AN_amr_female");
        bcf_info_t* AC_amr_female  = bcf_get_info(hdr, rec, "AC_amr_female");
        bcf_info_t* AF_amr_female  = bcf_get_info(hdr, rec, "AF_amr_female");
        int AC_amr_female2 = AC_amr_female->v1.i;
        int AN_amr_female2 = AN_amr_female->v1.i;
        float AF_amr_female2;

        if(AF_amr_female == NULL){
            if(AN_amr_female2 != 0){
                AF_amr_female2 = AC_amr_female2 / AN_amr_female2;
            }
            else{
                AF_amr_female2 = 0;
            }
        }else{
            AF_amr_female2 = AF_amr_female->v1.f;
        }

        out_file_stream << "AC_amr_female=" << AC_amr_female2 <<";";
        out_file_stream << "AN_amr_female=" << AN_amr_female2 <<";";
        out_file_stream << "AF_amr_female=" << AF_amr_female2 << ";";        

	//EAS_MALE-------------------------------------

	bcf_info_t* AN_eas_male  = bcf_get_info(hdr, rec, "AN_eas_male");
        bcf_info_t* AC_eas_male  = bcf_get_info(hdr, rec, "AC_eas_male");
        bcf_info_t* AF_eas_male  = bcf_get_info(hdr, rec, "AF_eas_male");
        int AC_eas_male2 = AC_eas_male->v1.i;
        int AN_eas_male2 = AN_eas_male->v1.i;
        float AF_eas_male2;

        if(AF_eas_male == NULL){
            if(AN_eas_male2 != 0){
                AF_eas_male2 = AC_eas_male2 / AN_eas_male2;
            }
            else{
                AF_eas_male2 = 0;
            }
        }else{
            AF_eas_male2 = AF_eas_male->v1.f;
        }

        out_file_stream << "AC_eas_male=" << AC_eas_male2 <<";";
        out_file_stream << "AN_eas_male=" << AN_eas_male2 <<";";
        out_file_stream << "AF_eas_male=" << AF_eas_male2 << ";";

	//EAS_FEMALE--------------------------------
	bcf_info_t* AN_eas_female  = bcf_get_info(hdr, rec, "AN_eas_female");
        bcf_info_t* AC_eas_female  = bcf_get_info(hdr, rec, "AC_eas_female");
        bcf_info_t* AF_eas_female  = bcf_get_info(hdr, rec, "AF_eas_female");
        int AC_eas_female2 = AC_eas_female->v1.i;
        int AN_eas_female2 = AN_eas_female->v1.i;
        float AF_eas_female2;

        if(AF_eas_female == NULL){
            if(AN_eas_female2 != 0){
                AF_eas_female2 = AC_eas_female2 / AN_eas_female2;
            }
            else{
                AF_eas_female2 = 0;
            }
        }else{
            AF_eas_female2 = AF_eas_female->v1.f;
        }

        out_file_stream << "AC_eas_female=" << AC_eas_female2 <<";";
        out_file_stream << "AN_eas_female=" << AN_eas_female2 <<";";
        out_file_stream << "AF_eas_female=" << AF_eas_female2 << ";";

	//NFE_MALE-------------------------------
	 bcf_info_t* AN_nfe_male  = bcf_get_info(hdr, rec, "AN_nfe_male");
        bcf_info_t* AC_nfe_male  = bcf_get_info(hdr, rec, "AC_nfe_male");
        bcf_info_t* AF_nfe_male  = bcf_get_info(hdr, rec, "AF_nfe_male");
        int AC_nfe_male2 = AC_nfe_male->v1.i;
        int AN_nfe_male2 = AN_nfe_male->v1.i;
        float AF_nfe_male2;

        if(AF_nfe_male == NULL){
            if(AN_nfe_male2 != 0){
                AF_nfe_male2 = AC_nfe_male2 / AN_nfe_male2;
            }
            else{
                AF_nfe_male2 = 0;
            }
        }else{
            AF_nfe_male2 = AF_nfe_male->v1.f;
        }

        out_file_stream << "AC_nfe_male=" << AC_nfe_male2 <<";";
        out_file_stream << "AN_nfe_male=" << AN_nfe_male2 <<";";
        out_file_stream << "AF_nfe_male=" << AF_nfe_male2 <<";";


	//NFE_FEMALE---------------------------
	bcf_info_t* AN_nfe_female  = bcf_get_info(hdr, rec, "AN_nfe_female");
        bcf_info_t* AC_nfe_female  = bcf_get_info(hdr, rec, "AC_nfe_female");
        bcf_info_t* AF_nfe_female  = bcf_get_info(hdr, rec, "AF_nfe_female");
        int AC_nfe_female2 = AC_nfe_female->v1.i;
        int AN_nfe_female2 = AN_nfe_female->v1.i;
        float AF_nfe_female2;

        if(AF_nfe_female == NULL){
            if(AN_nfe_female2 != 0){
                AF_nfe_female2 = AC_nfe_female2 / AN_nfe_female2;
            }
            else{
                AF_nfe_female2 = 0;
            }
        }else{
            AF_nfe_female2 = AF_nfe_female->v1.f;
        }

        out_file_stream << "AC_nfe_female=" << AC_nfe_female2 <<";";
        out_file_stream << "AN_nfe_female=" << AN_nfe_female2 <<";";
        out_file_stream << "AF_nfe_female=" << AF_nfe_female2 <<";";

	//ASJ_MALE-----------------------------
	bcf_info_t* AN_asj_male  = bcf_get_info(hdr, rec, "AN_asj_male");
        bcf_info_t* AC_asj_male  = bcf_get_info(hdr, rec, "AC_asj_male");
        bcf_info_t* AF_asj_male  = bcf_get_info(hdr, rec, "AF_asj_male");
        int AC_asj_male2 = AC_asj_male->v1.i;
        int AN_asj_male2 = AN_asj_male->v1.i;
        float AF_asj_male2;

        if(AF_asj_male == NULL){
            if(AN_asj_male2 != 0){
                AF_asj_male2 = AC_asj_male2 / AN_asj_male2;
            }
            else{
                AF_asj_male2 = 0;
            }
        }else{
            AF_asj_male2 = AF_asj_male->v1.f;
        }

        out_file_stream << "AC_asj_male=" << AC_asj_male2 <<";";
        out_file_stream << "AN_asj_male=" << AN_asj_male2 <<";";
        out_file_stream << "AF_asj_male=" << AF_asj_male2 <<";";


	//ASJ_FEMALE-----------------------------
	bcf_info_t* AN_asj_female  = bcf_get_info(hdr, rec, "AN_asj_female");
        bcf_info_t* AC_asj_female  = bcf_get_info(hdr, rec, "AC_asj_female");
        bcf_info_t* AF_asj_female  = bcf_get_info(hdr, rec, "AF_asj_female");
        int AC_asj_female2 = AC_asj_female->v1.i;
        int AN_asj_female2 = AN_asj_female->v1.i;
        float AF_asj_female2;

        if(AF_asj_female == NULL){
            if(AN_asj_female2 != 0){
                AF_asj_female2 = AC_asj_female2 / AN_asj_female2;
            }
            else{
                AF_asj_female2 = 0;
            }
        }else{
            AF_asj_female2 = AF_asj_female->v1.f;
        }

        out_file_stream << "AC_asj_female=" << AC_asj_female2 <<";";
        out_file_stream << "AN_asj_female=" << AN_asj_female2  <<";";
        out_file_stream << "AF_asj_female=" << AF_asj_female2  <<";";


	//OTH_MALE-------------------------------
	bcf_info_t* AN_oth_male  = bcf_get_info(hdr, rec, "AN_oth_male");
        bcf_info_t* AC_oth_male  = bcf_get_info(hdr, rec, "AC_oth_male");
        bcf_info_t* AF_oth_male  = bcf_get_info(hdr, rec, "AF_oth_male");
        int AC_oth_male2 = AC_oth_male->v1.i;
        int AN_oth_male2 = AN_oth_male->v1.i;
        float AF_oth_male2;

        if(AF_oth_male == NULL){
            if(AN_oth_male2 != 0){
                AF_oth_male2 = AC_oth_male2 / AN_oth_male2;
            }
            else{
                AF_oth_male2 = 0;
            }
        }else{
            AF_oth_male2 = AF_oth_male->v1.f;
        }

        out_file_stream << "AC_oth_male=" << AC_oth_male2 <<";";
        out_file_stream << "AN_oth_male=" << AN_oth_male2 <<";";
        out_file_stream << "AF_oth_male=" << AF_oth_male2 <<";";

	//OTH_FEMALE------------------------------
	bcf_info_t* AN_oth_female  = bcf_get_info(hdr, rec, "AN_oth_female");
        bcf_info_t* AC_oth_female  = bcf_get_info(hdr, rec, "AC_oth_female");
        bcf_info_t* AF_oth_female  = bcf_get_info(hdr, rec, "AF_oth_female");
        int AC_oth_female2 = AC_oth_female->v1.i;
        int AN_oth_female2 = AN_oth_female->v1.i;
        float AF_oth_female2;

        if(AF_oth_female == NULL){
            if(AN_oth_female2 != 0){
                AF_oth_female2 = AC_oth_female2 / AN_oth_female2;
            }
            else{
                AF_oth_female2 = 0;
            }
        }else{
            AF_oth_female2 = AF_oth_female->v1.f;
        }

        out_file_stream << "AC_oth_female=" << AC_oth_female2 <<";";
        out_file_stream << "AN_oth_female=" << AN_oth_female2 <<";";
        out_file_stream << "AF_oth_female=" << AF_oth_female2 <<";";

	//FIN_MALE--------------------------
	bcf_info_t* AN_fin_male  = bcf_get_info(hdr, rec, "AN_fin_male");
        bcf_info_t* AC_fin_male  = bcf_get_info(hdr, rec, "AC_fin_male");
        bcf_info_t* AF_fin_male  = bcf_get_info(hdr, rec, "AF_fin_male");
        int AC_fin_male2 = AC_fin_male->v1.i;
        int AN_fin_male2 = AN_fin_male->v1.i;
        float AF_fin_male2;

        if(AF_fin_male == NULL){
            if(AN_fin_male2 != 0){
                AF_fin_male2 = AC_fin_male2 / AN_fin_male2;
            }
            else{
                AF_fin_male2 = 0;
            }
        }else{
            AF_fin_male2 = AF_fin_male->v1.f;
        }

        out_file_stream << "AC_fin_male=" << AC_fin_male2 <<";";
        out_file_stream << "AN_fin_male=" << AN_fin_male2 <<";";
        out_file_stream << "AF_fin_male=" << AF_fin_male2 <<";";

	//FIN_FEMALE---------------------
	bcf_info_t* AN_fin_female  = bcf_get_info(hdr, rec, "AN_fin_female");
        bcf_info_t* AC_fin_female  = bcf_get_info(hdr, rec, "AC_fin_female");
        bcf_info_t* AF_fin_female  = bcf_get_info(hdr, rec, "AF_fin_female");
        int AC_fin_female2 = AC_fin_female->v1.i;
        int AN_fin_female2 = AN_fin_female->v1.i;
        float AF_fin_female2;

        if(AF_fin_female == NULL){
            if(AN_fin_female2 != 0){
                AF_fin_female2 = AC_fin_female2 / AN_fin_female2;
            }
            else{
                AF_fin_female2 = 0;
            }
        }else{
            AF_fin_female2 = AF_fin_female->v1.f;
        }

        out_file_stream << "AC_fin_female=" << AC_fin_female2 <<";";
        out_file_stream << "AN_fin_female=" << AN_fin_female2 <<";";
        out_file_stream << "AF_fin_female=" << AF_fin_female2 <<";";
	




	//AFR_MALE------------------------------
	bcf_info_t* AN_afr_male  = bcf_get_info(hdr, rec, "AN_afr_male");
        bcf_info_t* AC_afr_male  = bcf_get_info(hdr, rec, "AC_afr_male");
        bcf_info_t* AF_afr_male  = bcf_get_info(hdr, rec, "AF_afr_male");
        int AC_afr_male2 = AC_afr_male->v1.i;
        int AN_afr_male2 = AN_afr_male->v1.i;
        float AF_afr_male2;

        if(AF_afr_male == NULL){
            if(AN_afr_male2 != 0){
                AF_afr_male2 = AC_afr_male2 / AN_afr_male2;
            }
            else{
                AF_afr_male2 = 0;
            }
        }else{
            AF_afr_male2 = AF_afr_male->v1.f;
        }

        out_file_stream << "AC_afr_male=" << AC_afr_male2 <<";";
        out_file_stream << "AN_afr_male=" << AN_afr_male2 <<";";
        out_file_stream << "AF_afr_male=" << AF_afr_male2 <<";";

	//AFR_FEMALE----------------------------------

	bcf_info_t* AN_afr_female  = bcf_get_info(hdr, rec, "AN_afr_female");
        bcf_info_t* AC_afr_female  = bcf_get_info(hdr, rec, "AC_afr_female");
        bcf_info_t* AF_afr_female  = bcf_get_info(hdr, rec, "AF_afr_female");
        int AC_afr_female2 = AC_afr_female->v1.i;
        int AN_afr_female2 = AN_afr_female->v1.i;
        float AF_afr_female2;

        if(AF_afr_female == NULL){
            if(AN_afr_female2 != 0){
                AF_afr_female2 = AC_afr_female2 / AN_afr_female2;
            }
            else{
                AF_afr_female2 = 0;
            }
        }else{
            AF_afr_female2 = AF_afr_female->v1.f;
        }

        out_file_stream << "AC_afr_female=" << AC_afr_female2 <<";";
        out_file_stream << "AN_afr_female=" << AN_afr_female2 <<";";
        out_file_stream << "AF_afr_female=" << AF_afr_female2 <<";";


	//SAS_MALE----------------------------------
	bcf_info_t* AN_sas_male  = bcf_get_info(hdr, rec, "AN_sas_male");
        bcf_info_t* AC_sas_male  = bcf_get_info(hdr, rec, "AC_sas_male");
        bcf_info_t* AF_sas_male  = bcf_get_info(hdr, rec, "AF_sas_male");
        int AC_sas_male2 = AC_sas_male->v1.i;
        int AN_sas_male2 = AN_sas_male->v1.i;
        float AF_sas_male2;

        if(AF_sas_male == NULL){
            if(AN_sas_male2 != 0){
                AF_sas_male2 = AC_sas_male2 / AN_sas_male2;
            }
            else{
                AF_sas_male2 = 0;
            }
        }else{
            AF_sas_male2 = AF_sas_male->v1.f;
        }

        out_file_stream << "AC_sas_male=" << AC_sas_male2 <<";";
        out_file_stream << "AN_sas_male=" << AN_sas_male2 <<";";
        out_file_stream << "AF_sas_male=" << AF_sas_male2 <<";";

	//SAS_FEMALE---------------------------------
	bcf_info_t* AN_sas_female  = bcf_get_info(hdr, rec, "AN_sas_female");
        bcf_info_t* AC_sas_female  = bcf_get_info(hdr, rec, "AC_sas_female");
        bcf_info_t* AF_sas_female  = bcf_get_info(hdr, rec, "AF_sas_female");
        int AC_sas_female2 = AC_sas_female->v1.i;
        int AN_sas_female2 = AN_sas_female->v1.i;
        float AF_sas_female2;

        if(AF_sas_female == NULL){
            if(AN_sas_female2 != 0){
                AF_sas_female2 = AC_sas_female2 / AN_sas_female2;
            }
            else{
                AF_sas_female2 = 0;
            }
        }else{
            AF_sas_female2 = AF_sas_female->v1.f;
        }

        out_file_stream << "AC_sas_female=" << AC_sas_female2 <<";";
        out_file_stream << "AN_sas_female=" << AN_sas_female2 <<";";
        out_file_stream << "AF_sas_female=" << AF_sas_female2 <<";";

        //Southern European ancestry
        bcf_info_t* AN_nfe_seu  = bcf_get_info(hdr, rec, "AN_nfe_seu");
        bcf_info_t* AC_nfe_seu  = bcf_get_info(hdr, rec, "AC_nfe_seu");
        bcf_info_t* AF_nfe_seu  = bcf_get_info(hdr, rec, "AF_nfe_seu");
        int AC_nfe_seu2 = AC_nfe_seu->v1.i;
        int AN_nfe_seu2 = AN_nfe_seu->v1.i;
        float AF_nfe_seu2;

        if(AF_nfe_seu == NULL){
            if(AN_nfe_seu2 != 0){
                AF_nfe_seu2 = AC_nfe_seu2 / AN_nfe_seu2;
            }
            else{
                AF_nfe_seu2 = 0;
            }
        }else{
            AF_nfe_seu2 = AF_nfe_seu->v1.f;
        }

        out_file_stream << "AC_nfe_seu=" << AC_nfe_seu2 <<";";
        out_file_stream << "AN_nfe_seu=" << AN_nfe_seu2 <<";";
        out_file_stream << "AF_nfe_seu=" << AF_nfe_seu2 <<";";

        //Bulgarian (Eastern European) ancestry
        bcf_info_t* AN_nfe_bgr  = bcf_get_info(hdr, rec, "AN_nfe_bgr");
        bcf_info_t* AC_nfe_bgr  = bcf_get_info(hdr, rec, "AC_nfe_bgr");
        bcf_info_t* AF_nfe_bgr  = bcf_get_info(hdr, rec, "AF_nfe_bgr");
        int AC_nfe_bgr2 = AC_nfe_bgr->v1.i;
        int AN_nfe_bgr2 = AN_nfe_bgr->v1.i;
        float AF_nfe_bgr2;

        if(AF_nfe_bgr == NULL){
            if(AN_nfe_bgr2 != 0){
                AF_nfe_bgr2 = AC_nfe_bgr2 / AN_nfe_bgr2;
            }
            else{
                AF_nfe_bgr2 = 0;
            }
        }else{
            AF_nfe_bgr2 = AF_nfe_bgr->v1.f;
        }

        out_file_stream << "AC_nfe_bgr=" << AC_nfe_bgr2 <<";";
        out_file_stream << "AN_nfe_bgr=" << AN_nfe_bgr2 <<";";
        out_file_stream << "AF_nfe_bgr=" << AF_nfe_bgr2 <<";";


        //African-American/Afircan ancestry
        bcf_info_t* AN_afr  = bcf_get_info(hdr, rec, "AN_afr");
        bcf_info_t* AC_afr  = bcf_get_info(hdr, rec, "AC_afr");
        bcf_info_t* AF_afr  = bcf_get_info(hdr, rec, "AF_afr");
        int AC_afr2 = AC_afr->v1.i;
        int AN_afr2 = AN_afr->v1.i;
        float AF_afr2;

        if(AF_afr == NULL){
            if(AN_afr2 != 0){
                AF_afr2 = AC_afr2 / AN_afr2;
            }
            else{
                AF_afr2 = 0;
            }
        }else{
            AF_afr2 = AF_afr->v1.f;
        }

        out_file_stream << "AC_afr=" << AC_afr2 <<";";
        out_file_stream << "AN_afr=" << AN_afr2 <<";";
        out_file_stream << "AF_afr=" << AF_afr2 <<";";


        //South Asian ancestry
        bcf_info_t* AN_sas  = bcf_get_info(hdr, rec, "AN_sas");
        bcf_info_t* AC_sas  = bcf_get_info(hdr, rec, "AC_sas");
        bcf_info_t* AF_sas  = bcf_get_info(hdr, rec, "AF_sas");
        int AC_sas2 = AC_sas->v1.i;
        int AN_sas2 = AN_sas->v1.i;
        float AF_sas2;

        if(AF_sas == NULL){
            if(AN_sas2 != 0){
                AF_sas2 = AC_sas2 / AN_sas2;
            }
            else{
                AF_sas2 = 0;
            }
        }else{
            AF_sas2 = AF_sas->v1.f;
        }

        out_file_stream << "AC_sas=" << AC_sas2 <<";";
        out_file_stream << "AN_sas=" << AN_sas2 <<";";
        out_file_stream << "AF_sas=" << AF_sas2 <<";";


        //Other Non-Finnish European ancestry
        bcf_info_t* AN_nfe_onf  = bcf_get_info(hdr, rec, "AN_nfe_onf");
        bcf_info_t* AC_nfe_onf  = bcf_get_info(hdr, rec, "AC_nfe_onf");
        bcf_info_t* AF_nfe_onf  = bcf_get_info(hdr, rec, "AF_nfe_onf");
        int AC_nfe_onf2 = AC_nfe_onf->v1.i;
        int AN_nfe_onf2 = AN_nfe_onf->v1.i;
        float AF_nfe_onf2;

        if(AF_nfe_onf == NULL){
            if(AN_nfe_onf2 != 0){
                AF_nfe_onf2 = AC_nfe_onf2 / AN_nfe_onf2;
            }
            else{
                AF_nfe_onf2 = 0;
            }
        }else{
            AF_nfe_onf2 = AF_nfe_onf->v1.f;
        }

        out_file_stream << "AC_nfe_onf=" << AC_nfe_onf2 <<";";
        out_file_stream << "AN_nfe_onf=" << AN_nfe_onf2 <<";";
        out_file_stream << "AF_nfe_onf=" << AF_nfe_onf2 <<";";

        //Latino ancestry
        bcf_info_t* AN_amr  = bcf_get_info(hdr, rec, "AN_amr");
        bcf_info_t* AC_amr  = bcf_get_info(hdr, rec, "AC_amr");
        bcf_info_t* AF_amr  = bcf_get_info(hdr, rec, "AF_amr");
        int AC_amr2 = AC_amr->v1.i;
        int AN_amr2 = AN_amr->v1.i;
        float AF_amr2;

        if(AF_amr == NULL){
            if(AN_amr2 != 0){
                AF_amr2 = AC_amr2 / AN_amr2;
            }
            else{
                AF_amr2 = 0;
            }
        }else{
            AF_amr2 = AF_amr->v1.f;
        }

        out_file_stream << "AC_amr=" << AC_amr2 <<";";
        out_file_stream << "AN_amr=" << AN_amr2 <<";";
        out_file_stream << "AF_amr=" << AF_amr2 <<";";

        //East Asian ancestry
        bcf_info_t* AN_eas  = bcf_get_info(hdr, rec, "AN_eas");
        bcf_info_t* AC_eas  = bcf_get_info(hdr, rec, "AC_eas");
        bcf_info_t* AF_eas  = bcf_get_info(hdr, rec, "AF_eas");
        int AC_eas2 = AC_eas->v1.i;
        int AN_eas2 = AN_eas->v1.i;
        float AF_eas2;

        if(AF_eas == NULL){
            if(AN_eas2 != 0){
                AF_eas2 = AC_eas2 / AN_eas2;
            }
            else{
                AF_eas2 = 0;
            }
        }else{
            AF_eas2 = AF_eas->v1.f;
        }

        out_file_stream << "AC_eas=" << AC_eas2 <<";";
        out_file_stream << "AN_eas=" << AN_eas2 <<";";
        out_file_stream << "AF_eas=" << AF_eas2 <<";";


        //Swedish ancestry
        bcf_info_t* AN_nfe_swe  = bcf_get_info(hdr, rec, "AN_nfe_swe");
        bcf_info_t* AC_nfe_swe  = bcf_get_info(hdr, rec, "AC_nfe_swe");
        bcf_info_t* AF_nfe_swe  = bcf_get_info(hdr, rec, "AF_nfe_swe");
        int AC_nfe_swe2 = AC_nfe_swe->v1.i;
        int AN_nfe_swe2 = AN_nfe_swe->v1.i;
        float AF_nfe_swe2;

        if(AF_nfe_swe == NULL){
            if(AN_nfe_swe2 != 0){
                AF_nfe_swe2 = AC_nfe_swe2 / AN_nfe_swe2;
            }
            else{
                AF_nfe_swe2 = 0;
            }
        }else{
            AF_nfe_swe2 = AF_nfe_swe->v1.f;
        }

        out_file_stream << "AC_nfe_swe=" << AC_nfe_swe2 <<";";
        out_file_stream << "AN_nfe_swe=" << AN_nfe_swe2 <<";";
        out_file_stream << "AF_nfe_swe=" << AF_nfe_swe2 <<";";

        //North-Western European ancestry
        bcf_info_t* AN_nfe_nwe  = bcf_get_info(hdr, rec, "AN_nfe_nwe");
        bcf_info_t* AC_nfe_nwe  = bcf_get_info(hdr, rec, "AC_nfe_nwe");
        bcf_info_t* AF_nfe_nwe  = bcf_get_info(hdr, rec, "AF_nfe_nwe");
        int AC_nfe_nwe2 = AC_nfe_nwe->v1.i;
        int AN_nfe_nwe2 = AN_nfe_nwe->v1.i;
        float AF_nfe_nwe2;

        if(AF_nfe_nwe == NULL){
            if(AN_nfe_nwe2 != 0){
                AF_nfe_nwe2 = AC_nfe_nwe2 / AN_nfe_nwe2;
            }
            else{
                AF_nfe_nwe2 = 0;
            }
        }else{
            AF_nfe_nwe2 = AF_nfe_nwe->v1.f;
        }

        out_file_stream << "AC_nfe_nwe=" << AC_nfe_nwe2 <<";";
        out_file_stream << "AN_nfe_nwe=" << AN_nfe_nwe2 <<";";
        out_file_stream << "AF_nfe_nwe=" << AF_nfe_nwe2 <<";";


        //Japanese ancestry
        bcf_info_t* AN_eas_jpn  = bcf_get_info(hdr, rec, "AN_eas_jpn");
        bcf_info_t* AC_eas_jpn  = bcf_get_info(hdr, rec, "AC_eas_jpn");
        bcf_info_t* AF_eas_jpn  = bcf_get_info(hdr, rec, "AF_eas_jpn");
        int AC_eas_jpn2 = AC_eas_jpn->v1.i;
        int AN_eas_jpn2 = AN_eas_jpn->v1.i;
        float AF_eas_jpn2;

        if(AF_eas_jpn == NULL){
            if(AN_eas_jpn2 != 0){
                AF_eas_jpn2 = AC_eas_jpn2 / AN_eas_jpn2;
            }
            else{
                AF_eas_jpn2 = 0;
            }
        }else{
            AF_eas_jpn2 = AF_eas_jpn->v1.f;
        }

        out_file_stream << "AC_eas_jpn=" << AC_eas_jpn2 <<";";
        out_file_stream << "AN_eas_jpn=" << AN_eas_jpn2 <<";";
        out_file_stream << "AF_eas_jpn=" << AF_eas_jpn2 <<";";

        //Korean ancestry
        bcf_info_t* AN_eas_kor  = bcf_get_info(hdr, rec, "AN_eas_kor");
        bcf_info_t* AC_eas_kor  = bcf_get_info(hdr, rec, "AC_eas_kor");
        bcf_info_t* AF_eas_kor  = bcf_get_info(hdr, rec, "AF_eas_kor");
        int AC_eas_kor2 = AC_eas_kor->v1.i;
        int AN_eas_kor2 = AN_eas_kor->v1.i;
        float AF_eas_kor2;

        if(AF_eas_kor == NULL){
            if(AN_eas_kor2 != 0){
                AF_eas_kor2 = AC_eas_kor2 / AN_eas_kor2;
            }
            else{
                AF_eas_kor2 = 0;
            }
        }else{
            AF_eas_kor2 = AF_eas_kor->v1.f;
        }

        out_file_stream << "AC_eas_kor=" << AC_eas_kor2 <<";";
        out_file_stream << "AN_eas_kor=" << AN_eas_kor2 <<";";
        out_file_stream << "AF_eas_kor=" << AF_eas_kor2 <<";";

        //Other East Asian ancestry
        bcf_info_t* AN_eas_oea  = bcf_get_info(hdr, rec, "AN_eas_oea");
        bcf_info_t* AC_eas_oea  = bcf_get_info(hdr, rec, "AC_eas_oea");
        bcf_info_t* AF_eas_oea  = bcf_get_info(hdr, rec, "AF_eas_oea");
        int AC_eas_oea2 = AC_eas_oea->v1.i;
        int AN_eas_oea2 = AN_eas_oea->v1.i;
        float AF_eas_oea2;

        if(AF_eas_oea == NULL){
            if(AN_eas_oea2 != 0){
                AF_eas_oea2 = AC_eas_oea2 / AN_eas_oea2;
            }
            else{
                AF_eas_oea2 = 0;
            }
        }else{
            AF_eas_oea2 = AF_eas_oea->v1.f;
        }

        out_file_stream << "AC_eas_oea=" << AC_eas_oea2 <<";";
        out_file_stream << "AN_eas_oea=" << AN_eas_oea2 <<";";
        out_file_stream << "AF_eas_oea=" << AF_eas_oea2 <<";";

        //Estonian ancestry
        bcf_info_t* AN_nfe_est  = bcf_get_info(hdr, rec, "AN_nfe_est");
        bcf_info_t* AC_nfe_est  = bcf_get_info(hdr, rec, "AC_nfe_est");
        bcf_info_t* AF_nfe_est  = bcf_get_info(hdr, rec, "AF_nfe_est");
        int AC_nfe_est2 = AC_nfe_est->v1.i;
        int AN_nfe_est2 = AN_nfe_est->v1.i;
        float AF_nfe_est2;

        if(AF_nfe_est == NULL){
            if(AN_nfe_est2 != 0){
                AF_nfe_est2 = AC_nfe_est2 / AN_nfe_est2;
            }
            else{
                AF_nfe_est2 = 0;
            }
        }else{
            AF_nfe_est2 = AF_nfe_est->v1.f;
        }

        out_file_stream << "AC_nfe_est=" << AC_nfe_est2 <<";";
        out_file_stream << "AN_nfe_est=" << AN_nfe_est2 <<";";
        out_file_stream << "AF_nfe_est=" << AF_nfe_est2 <<";";


        //Non-Finnish European ancestry
        bcf_info_t* AN_nfe  = bcf_get_info(hdr, rec, "AN_nfe");
        bcf_info_t* AC_nfe  = bcf_get_info(hdr, rec, "AC_nfe");
        bcf_info_t* AF_nfe  = bcf_get_info(hdr, rec, "AF_nfe");
        int AC_nfe2 = AC_nfe->v1.i;
        int AN_nfe2 = AN_nfe->v1.i;
        float AF_nfe2;

        if(AF_nfe == NULL){
            if(AN_nfe2 != 0){
                AF_nfe2 = AC_nfe2 / AN_nfe2;
            }
            else{
                AF_nfe2 = 0;
            }
        }else{
            AF_nfe2 = AF_nfe->v1.f;
        }

        out_file_stream << "AC_nfe=" << AC_nfe2 <<";";
        out_file_stream << "AN_nfe=" << AN_nfe2 <<";";
        out_file_stream << "AF_nfe=" << AF_nfe2 <<";";

        //Finnish ancestry
        bcf_info_t* AN_fin  = bcf_get_info(hdr, rec, "AN_fin");
        bcf_info_t* AC_fin  = bcf_get_info(hdr, rec, "AC_fin");
        bcf_info_t* AF_fin  = bcf_get_info(hdr, rec, "AF_fin");
        int AC_fin2 = AC_fin->v1.i;
        int AN_fin2 = AN_fin->v1.i;
        float AF_fin2;

        if(AF_fin == NULL){
            if(AN_fin2 != 0){
                AF_fin2 = AC_fin2 / AN_fin2;
            }
            else{
                AF_fin2 = 0;
            }
        }else{
            AF_fin2 = AF_fin->v1.f;
        }

        out_file_stream << "AC_fin=" << AC_fin2 <<";";
        out_file_stream << "AN_fin=" << AN_fin2 <<";";
        out_file_stream << "AF_fin=" << AF_fin2 <<";";

        //Ashkenazi Jewish ancestry
        bcf_info_t* AN_asj  = bcf_get_info(hdr, rec, "AN_asj");
        bcf_info_t* AC_asj  = bcf_get_info(hdr, rec, "AC_asj");
        bcf_info_t* AF_asj  = bcf_get_info(hdr, rec, "AF_asj");
        int AC_asj2 = AC_asj->v1.i;
        int AN_asj2 = AN_asj->v1.i;
        float AF_asj2;

        if(AF_asj == NULL){
            if(AN_asj2 != 0){
                AF_asj2 = AC_asj2 / AN_asj2;
            }
            else{
                AF_asj2 = 0;
            }
        }else{
            AF_asj2 = AF_asj->v1.f;
        }

        out_file_stream << "AC_asj=" << AC_asj2 <<";";
        out_file_stream << "AN_asj=" << AN_asj2 <<";";
        out_file_stream << "AF_asj=" << AF_asj2 <<";";

        //Other ancestry
        bcf_info_t* AN_oth  = bcf_get_info(hdr, rec, "AN_oth");
        bcf_info_t* AC_oth  = bcf_get_info(hdr, rec, "AC_oth");
        bcf_info_t* AF_oth  = bcf_get_info(hdr, rec, "AF_oth");
        int AC_oth2 = AC_oth->v1.i;
        int AN_oth2 = AN_oth->v1.i;
        float AF_oth2;

        if(AF_oth == NULL){
            if(AN_oth2 != 0){
                AF_oth2 = AC_oth2 / AN_oth2;
            }
            else{
                AF_oth2 = 0;
            }
        }else{
            AF_oth2 = AF_oth->v1.f;
        }

        out_file_stream << "AC_oth=" << AC_oth2 <<";";
        out_file_stream << "AN_oth=" << AN_oth2 <<";";
        out_file_stream << "AF_oth=" << AF_oth2 <<";";

	out_file_stream << vep_val;    
	out_file_stream << std::endl;
    }



    bcf_hdr_destroy(hdr);
    bcf_close(inf);
    bcf_destroy(rec);
}













