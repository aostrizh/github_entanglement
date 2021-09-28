#include <iostream>
#include <fstream>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>

#include "like_ivashkin_wants_it.h"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>

#include <TChain.h>
#include <TObjArray.h>
#include <TString.h>
#include <TDirectory.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>
#include <TGraphErrors.h>
#include <TObjString.h>
#include "../calo_analysis/Readme_reader.h"
#include "ChannelEntry.h"
using namespace std;

#define UseManyRoots 0
#define PresentDiffuser 1
#define CalculateRatio 1
#define DrawTime 1
#define UseIntegralCut 1
#define UseNotEntangledPhotons 0
#define TotalEnergyCut 0
#define DrawToPDF 0
#define calculate_CHSH 1
/////////////////////////////////////////////////////
#if calculate_CHSH
/////////////////Calculate global CHSH (global N_i Summ of different polarizers N_i)
void count_coincidences(Int_t n, Int_t NumEvents[16][16], Int_t &N_parallel, Int_t &N_perpendicular, Int_t &N_mixed_1, Int_t &N_mixed_2, TString string)
{
        N_parallel=0;
        N_perpendicular = 0;
        N_mixed_1 = 0;
        N_mixed_2 = 0;  

    if (string == "clockwise" || string == "counterclockwise")
    {
        for (Int_t a = 0; a < 16; a++)
        {   
            Int_t b = 0;
            if (string == "clockwise") b = a+n;
            if (string == "counterclockwise") b = a-n+16;
            Int_t a_1 = a+4;
            Int_t b_1 = b+4;
            while (a > 15) a=a-16;
            while (b > 15) b=b-16;

            while (a_1> 15) a_1=a_1-16;
            while (b_1 > 15) b_1=b_1-16;    
            N_parallel += NumEvents[a][b];
            N_perpendicular += NumEvents[a_1][b_1];
            N_mixed_1 += NumEvents[a][b_1];
            N_mixed_2 += NumEvents[a_1][b];        
        }        
    }
    if (string == "both") 
    {
         for (Int_t a = 0; a < 16; a++)
        {
            Int_t a_1 = a+4;

            Int_t b = a+n;
            Int_t b_1 = b+4;

            Int_t b_ccl = a-n+16;
            Int_t b_ccl_1 = b_ccl +4;

            while (a > 15) a=a-16;
            while (a_1> 15) a_1=a_1-16;

            while (b > 15) b=b-16;
            while (b_1 > 15) b_1=b_1-16;  

            while (b_ccl > 15) b_ccl=b_ccl-16;
            while (b_ccl_1 > 15) b_ccl_1=b_ccl_1-16;
  
            N_parallel += NumEvents[a][b]; N_parallel += NumEvents[a][b_ccl];
            N_perpendicular += NumEvents[a_1][b_1];  N_perpendicular += NumEvents[a_1][b_ccl_1];
            N_mixed_1 += NumEvents[a][b_1]; N_mixed_1 += NumEvents[a][b_ccl_1];
            N_mixed_2 += NumEvents[a_1][b]; N_mixed_2 += NumEvents[a_1][b_ccl]; 
        }
    }

}

Float_t global_E_coeff(Int_t n, Int_t NumEvents[16][16], TString string)
{
    Int_t N_parallel = 0; Int_t N_perpendicular = 0;Int_t N_mixed_1 = 0;Int_t N_mixed_2 = 0;
    count_coincidences(n,NumEvents, N_parallel, N_perpendicular, N_mixed_1, N_mixed_2, string);
    return (float)(N_parallel+N_perpendicular-N_mixed_1-N_mixed_2)/(N_parallel+N_perpendicular+N_mixed_1+N_mixed_2)
    ;
}
Float_t global_sqr_E_error(Int_t n,Int_t NumEvents[16][16], TString string) //calculating square error of correlation coefficient
{
    Int_t N_parallel = 0; Int_t N_perpendicular = 0;Int_t N_mixed_1 = 0;Int_t N_mixed_2 = 0;
    count_coincidences(n,NumEvents, N_parallel, N_perpendicular, N_mixed_1, N_mixed_2, string);
    return (float)4.0/pow((N_parallel+N_perpendicular+N_mixed_1+N_mixed_2),4)*(pow(N_mixed_1+N_mixed_2,2)*(N_parallel+N_perpendicular)+pow(N_parallel+N_perpendicular,2)*(N_mixed_1+N_mixed_2));
}

Float_t global_calculate_CHSH(Int_t n,Int_t NumEvents[16][16], TString string)
{
    return
    3*global_E_coeff(n, NumEvents, string) - global_E_coeff(3*n,NumEvents, string)
    ;
}

Float_t global_calculate_CHSH_error(Int_t n,Int_t NumEvents[16][16], TString string)
{
    return
    global_sqr_E_error(n, NumEvents, string)*3 + global_sqr_E_error(3*n, NumEvents, string)
    ;
}

///////////////////////
///////////////////////Calculating average CHSH for all polarizers 

Float_t E_coeff(Int_t NumEvents[16][16], Int_t a, Int_t n, TString string)
{
    Int_t b = 0;
    if (string == "clockwise") 
    b = a+n;
    if (string == "counterclockwise") b = a-n+16;

    Int_t a_1 = a+4;
    Int_t b_1 = b+4;
    while (a > 15) a=a-16;
    while (b > 15) b=b-16;

    while (a_1> 15) a_1=a_1-16;
    while (b_1 > 15) b_1=b_1-16;    
    Int_t N_parallel=NumEvents[a][b];
    Int_t N_perpendicular = NumEvents[a_1][b_1];
    Int_t N_mixed_1 = NumEvents[a][b_1];
    Int_t N_mixed_2 = NumEvents[a_1][b];
    return (float)(N_parallel+N_perpendicular-N_mixed_1-N_mixed_2)/(N_parallel+N_perpendicular+N_mixed_1+N_mixed_2)
    ;
}
Float_t sqr_E_error(Int_t NumEvents[16][16], Int_t a, Int_t n, TString string) //calculating square error of correlation coefficient
{
    Int_t b = 0;
    if (string == "clockwise") 
    b = a+n;
    if (string == "counterclockwise") b = a-n+16;
    Int_t a_1 = a+4;
    Int_t b_1 = b+4;
    while (a > 15) a=a-16;
    while (b > 15) b=b-16;

    while (a_1> 15) a_1=a_1-16;
    while (b_1 > 15) b_1=b_1-16;    
    Int_t N_parallel=NumEvents[a][b];
    Int_t N_perpendicular = NumEvents[a_1][b_1];
    Int_t N_mixed_1 = NumEvents[a][b_1];
    Int_t N_mixed_2 = NumEvents[a_1][b];
    //return pow(2*(N_mixed_1+N_mixed_2)/pow((N_parallel+N_perpendicular+N_mixed_1+N_mix)ed_2),2),2)*(N_parallel+N_perpendicular+N_mixed_1+N_mixed_2);
    return 
    (float)4.0/pow((N_parallel+N_perpendicular+N_mixed_1+N_mixed_2),4)*(pow(N_mixed_1+N_mixed_2,2)*(N_parallel+N_perpendicular)+pow(N_parallel+N_perpendicular,2)*(N_mixed_1+N_mixed_2))

;
}


Float_t calculate_local_CHSH(Int_t NumEvents[16][16], Int_t a_angle, Int_t angle_between_counters, TString string)
{
    return
    E_coeff(NumEvents,a_angle,angle_between_counters, string) -
    E_coeff(NumEvents,a_angle,angle_between_counters, string) +
    E_coeff(NumEvents,a_angle+2*angle_between_counters,angle_between_counters, string) +
    E_coeff(NumEvents,a_angle+2*angle_between_counters,3*angle_between_counters, string);
}

Float_t calculate_local_CHSH_error(Int_t NumEvents[16][16], Int_t a_angle, Int_t angle_between_counters, TString string)
{

    return
    sqr_E_error(NumEvents,a_angle,angle_between_counters, string) +
    sqr_E_error(NumEvents,a_angle,3*angle_between_counters, string) +
    sqr_E_error(NumEvents,a_angle+2*angle_between_counters,angle_between_counters, string) +
    sqr_E_error(NumEvents,a_angle+2*angle_between_counters,3*angle_between_counters, string);

    
}



#endif

/////////////////////////////////////////////////////
void ________FINAL_ANALYSIS_WRONG_WAY_______()
{
	gStyle->SetOptFit(1);
	//gStyle->SetOptStat(1111);
	TString source_path = "/home/doc/entanglement/root_files_data/with_scatterer/big_file/";
#if UseNotEntangledPhotons
    TString result_path  = source_path + "result_ev_b_ev_decoh_new";
#else
    TString result_path = source_path + "result_ev_b_ev_entangled_new";
#endif
    TFile *result_root = new TFile (result_path+".root", "RECREATE");

    const Int_t left_NaI_range = 180;
    const Int_t righr_NaI_range = 280;
    Float_t interval_width = 1.35;

////////////////
    Int_t low_cut[32] = {0}; Int_t high_cut[32] = {0}; Float_t total_low_cut[32] = {0}; Float_t total_high_cut[32] = {0};
    Float_t sigma_i[32] = {0}; Float_t mean_int[32] = {0};

    Int_t low_amp_cut[32] = {0};
    Int_t high_amp_cut[32] = {0};

    Float_t average_scatterer_peak_position[32] = {0};
    Float_t sigma_scatterer[32] = {0};
    Float_t low_energy_cut[32] = {0};
    Float_t high_energy_cut[32] = {0};

    Float_t sum_energy[32] = {0};
    Float_t sigma_energy[32] = {0};
    Float_t low_energy_sum[32] = {0};
    Float_t high_energy_sum[32] = {0};

/////////////////////
        const Int_t histos_bins_number = 200;
        const Int_t histos_left_boarder = 1;
        const Int_t histos_right_boarder = 450;
        TH1F *true_sum_hist_left = new TH1F ("TRUE_peak_left_sum_spectrum","TRUE_peak_left_sum_spectrum",histos_bins_number,histos_left_boarder,histos_right_boarder);
        TH1F *true_sum_hist_right = new TH1F ("TRUE_peak_right_sum_spectrum","TRUE_peak_right_sum_spectrum",histos_bins_number,histos_left_boarder,histos_right_boarder);
        TH1F *true_sc_sum_hist_left = new TH1F ("TRUE_peak_scatterer_left_sum_spectrum","TRUE_peak_scatterer_left_sum_spectrum",200,1,700);
        TH1F *true_sc_sum_hist_right = new TH1F ("TRUE_peak_scatterer_right_sum_spectrum","TRUE_peak_scatterer_right_sum_spectrum",200,1,700);
/////////////////////
    Int_t channel = 0;
    Int_t Board_Ch = 0;
	Int_t global_counter = 1;

    Readme_reader *readme_reader_ptr = new Readme_reader("");
	readme_reader_ptr->SetInfoFile(source_path + "/" + "run_info.txt");

    TH1F *hist_charge[16][16];

    TH1F *charge_hist[32];
    Int_t NumEvents[16][16] = {0};     
    Float_t angle_arr[16];
    Float_t angle_arr_err[16];


    for (Int_t channel_number = 0; channel_number < 16; channel_number++)
    {
        for (Int_t channel_number_2 = 16; channel_number_2 < 32; channel_number_2++)    
        {
            hist_charge[channel_number][channel_number_2-16] = 
            new TH1F(Form("integral_in_gate_%i_%i",channel_number,channel_number_2),
            Form("integral_in_gate_%i_%i",channel_number,channel_number_2),140,10,1000);   
        }     
    }    

    for (Int_t i = 0; i <16; i++)
    {
        angle_arr[i] = (float)i*22.5;

        angle_arr_err[i] = 0.;
    }
    Float_t total_events_for_angle_diff[16] = {0};
    Float_t total_events_for_angle_diff_err[16] = {0};

    TChain *PMT_tree = new TChain;

///////////////////////////////////////////////
#if UseManyRoots
	TObjArray files_names;
	TSystemDirectory dir(source_path, source_path);
	TList *files = dir.GetListOfFiles();
	if(files) 
	{
		TSystemFile *file; 
		TString fname; 
		TIter next(files); 
		while((file=(TSystemFile*)next())) 
		{ 
			fname = file->GetName(); 
			if(!file->IsDirectory() && fname.EndsWith("07a8de9a_20210224_164743.root")) 
			{
				files_names.Add(new TObjString(fname));
				PMT_tree->AddFile( (source_path + fname + "/adc64_data").Data() );
			}
		}
	}
    Int_t total_files = files_names.GetLast()+1;
	printf("Total files: %i\n", total_files);
//////////////////
#else
		PMT_tree->AddFile( (source_path + "calibrated_time.root" + "/adc64_data").Data() );

#endif

///////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////
///////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////
///////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\////////////////////////
    const Int_t calculate_Events_Number =  
    
    PMT_tree->GetEntries()/100
    ;
    const Int_t events_wanted = 10000000;
    Int_t Events_for_cuts = 0;
#if UseIntegralCut    
    Events_for_cuts = events_wanted;
    if (calculate_Events_Number < events_wanted) Events_for_cuts = calculate_Events_Number;
    if ((int)calculate_Events_Number/2 > 20000000) Events_for_cuts = (int)calculate_Events_Number/1;
#endif
///////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////
///////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////
///////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\////////////////////////



#if PresentDiffuser    
        short_energy_ChannelEntry *short_channel_info = new short_energy_ChannelEntry[35];
        for(Int_t ch = 0; ch < 35; ch++)
	        (short_channel_info+ch)->SetBranch(PMT_tree, ch);
#else
        short_energy_ChannelEntry *short_channel_info = new short_energy_ChannelEntry[34];
	    for(Int_t ch = 0; ch < 34; ch++)
	        (short_channel_info+ch)->SetBranch(PMT_tree, ch);
#endif
/////////////////////////////
    /////////////////////////integral_in_gate_cut_range_cuts

        TH1F *peak_histo[35];
        TH1F *peak_histo_w_o_cuts[35];
        
        for (Int_t channel_number = 0; channel_number < 32; channel_number++)
        {
            TString hernya = Form("integral_in_gate_channel_%i",channel_number);
            peak_histo[channel_number] = new TH1F (hernya.Data(),hernya.Data(),histos_bins_number,histos_left_boarder,histos_right_boarder);
        }
        for (Int_t iEvent = 0; iEvent < Events_for_cuts; iEvent++)
        {
            PMT_tree->GetEntry(iEvent);
            for (Int_t channel_number = 0; channel_number < 32; channel_number++)
            {
                Int_t sc_number = 33;
                if (channel_number < 16) sc_number = 32;
                if (


                abs(short_channel_info[33].peak_pos - short_channel_info[32].peak_pos) < 30

                &&   short_channel_info[channel_number].peak_pos > 1600
                &&   short_channel_info[33].peak_pos > 1600
                &&   short_channel_info[32].peak_pos > 1600
                && (short_channel_info[channel_number].peak_pos - short_channel_info[sc_number].peak_pos) > 100
                && (short_channel_info[channel_number].peak_pos - short_channel_info[sc_number].peak_pos) < 200


                && short_channel_info[33].amp < 60000
                && short_channel_info[33].amp > 400
                && short_channel_info[32].amp < 60000     
                && short_channel_info[32].amp > 400 
                && short_channel_info[channel_number].amp > 400 && short_channel_info[channel_number].amp < 60000
                && short_channel_info[channel_number].integral_in_gate > 4


#if UseIntegralCut

                && short_channel_info[32].integral_in_gate > 140
                && short_channel_info[32].integral_in_gate < 400             
#endif         
#if UseNotEntangledPhotons

                && short_channel_info[34].amp < 50000
                && short_channel_info[34].amp > 400
                && (short_channel_info[34].peak_pos - short_channel_info[32].peak_pos) < 60
                && (short_channel_info[34].peak_pos - short_channel_info[32].peak_pos) >-30                 
#else 
                && short_channel_info[34].amp < 400
                && abs(short_channel_info[34].peak_pos - short_channel_info[32].peak_pos) > 100
#endif

                ) 
                    {
                        peak_histo[channel_number]->Fill(short_channel_info[channel_number].integral_in_gate);
                    }

            }
        }

            for (Int_t channel_number = 0; channel_number < 32; channel_number++)
            {
				TF1 *gaus_fit = new TF1 ("gaus_func","gaus",left_NaI_range,righr_NaI_range);
				peak_histo[channel_number]->Fit("gaus_func","R");

				TF1 *gaus_2_fit = new TF1 ("gaus_func_2","gaus",
				(gaus_fit->GetParameter(1))-gaus_fit->GetParameter(2),
				(gaus_fit->GetParameter(1))+gaus_fit->GetParameter(2));
				peak_histo[channel_number]->Fit("gaus_func_2","R");

                mean_int[channel_number] = gaus_2_fit->GetParameter(1);
                sigma_i[channel_number] = gaus_2_fit->GetParameter(2);
				low_cut[channel_number] = gaus_2_fit->GetParameter(1)-interval_width*gaus_2_fit->GetParameter(2)
                
                ;
				high_cut[channel_number] = gaus_2_fit->GetParameter(1)+interval_width*gaus_2_fit->GetParameter(2)
                
                ;
                hist_like_ivashkin_wants_it(peak_histo[channel_number],"Energy [keV]", "Counts");
                peak_histo[channel_number]->Write();

                if (channel_number < 16) true_sum_hist_left->Add(peak_histo[channel_number]);
                if (channel_number >=16) true_sum_hist_right->Add(peak_histo[channel_number]);
                #if DrawToPDF
                TCanvas *canv = new TCanvas("canv","canv");
                canv->cd();
                peak_histo[channel_number]->Draw();

                canv->SaveAs(result_path+".pdf(",".pdf");
                #endif

                delete peak_histo[channel_number];

            }
            TH1F *int_without_cuts = new TH1F("ch_34.integral_in_gate","ch_34.integral_in_gate",200,1,300);
            TH1F *int_with_cuts = new TH1F("ch_34.integral_in_gate_with_time_cuts","ch_34.integral_in_gate_with_time_cuts",200,1,300);

                TCanvas *canv_0 = new TCanvas("canv_0","canv_0");
                canv_0->Divide(2);
                canv_0->cd(1);
                true_sum_hist_left->Draw();
                canv_0->cd(2);
                true_sum_hist_right->Draw();


                canv_0->SaveAs(result_path+".pdf(",".pdf");
            
/*            #if DrawToPDF
            TCanvas *temp_canv = new TCanvas("temp_canv","temp_canv");
            temp_canv->Divide(2);
            temp_canv->cd(1);
            PMT_tree->Draw("channel_34.integral_in_gate >> ch_34.integral_in_gate","channel_34.amp > 400 && channel_34.amp < 60000");
            int_without_cuts->GetYaxis()->SetRangeUser(0,(Int_t)1.05*int_without_cuts->GetMaximum());
            temp_canv->cd(2);
            PMT_tree->Draw("channel_34.integral_in_gate >> ch_34.integral_in_gate_with_time_cuts","channel_34.amp > 400 && channel_34.amp < 60000 && channel_34.peak_pos - channel_32.peak_pos > -30 && channel_34.peak_pos - channel_32.peak_pos < 60");
            int_with_cuts->GetYaxis()->SetRangeUser(0,(Int_t)1.05*int_without_cuts->GetMaximum());

            temp_canv->SaveAs(result_path + ".pdf(",".pdf");
            #endif
*/            /////////////////////////////////////////////////


////////////////////////energy in diffuser
        for (Int_t channel_number = 0; channel_number < 32; channel_number++)
        {
            TString hernya = Form("integral_in_gate_diffuser_ch_cut_%i",channel_number);
            peak_histo[channel_number] = new TH1F (hernya.Data(),hernya.Data(),100,0,200);
        }


        for (Int_t iEvent = 0; iEvent < Events_for_cuts; iEvent++)
        {
            PMT_tree->GetEntry(iEvent);
            for (Int_t channel_number = 0; channel_number < 16; channel_number++)
            {

                if (short_channel_info[channel_number].integral_in_gate > low_cut[channel_number]
                && short_channel_info[channel_number].integral_in_gate < high_cut[channel_number]

                && (short_channel_info[channel_number].peak_pos - short_channel_info[32].peak_pos) > 100
                && (short_channel_info[channel_number].peak_pos - short_channel_info[32].peak_pos) < 200

                &&   short_channel_info[channel_number].peak_pos > 1600
                &&   short_channel_info[33].peak_pos > 1600
                &&   short_channel_info[32].peak_pos > 1600



                && short_channel_info[channel_number].amp > 400 && short_channel_info[channel_number].amp < 60000         

                && short_channel_info[33].amp < 60000
                && short_channel_info[33].amp > 400
                && short_channel_info[32].amp < 60000     
                && short_channel_info[32].amp > 400 
                && short_channel_info[channel_number].amp > 400 && short_channel_info[channel_number].amp < 60000


#if UseNotEntangledPhotons

                && short_channel_info[34].amp < 50000
                && short_channel_info[34].amp > 400
                && (short_channel_info[34].peak_pos - short_channel_info[32].peak_pos) < 60
                && (short_channel_info[34].peak_pos - short_channel_info[32].peak_pos) >-30 

#else 
                && short_channel_info[34].amp < 400
                && abs(short_channel_info[34].peak_pos - short_channel_info[32].peak_pos) > 100
#endif
      



                ) 
                    peak_histo[channel_number]->Fill(short_channel_info[34].integral_in_gate);
            }
        }

    TString hist_name = "integral_in_gate_diffuser";
    TH1F *d_hist = new TH1F (hist_name.Data(),hist_name.Data(),100,0,200);


            for (Int_t channel_number = 0; channel_number < 16; channel_number++)
            {


                d_hist->Add(peak_histo[channel_number]);

                delete peak_histo[channel_number];

            }
                hist_like_ivashkin_wants_it(d_hist,"Energy [keV]", "Counts");
                d_hist->Write();

                TCanvas *canv_1 = new TCanvas("canv","canv");
                canv_1-> cd();
                d_hist->Draw();
                canv_1->SaveAs(result_path+".pdf(",".pdf");




                



///////////////////////////////
/////////energy in scatterer
/////////full_spectrum
/*            TCanvas *canv_00 = new TCanvas("canv_00","canv_00");
            canv_00->Divide(2);
            canv_00->cd(1);
            PMT_tree->Draw("channel_32.integral_in_gate >> TRUE_peak_scatterer_left_sum_spectrum"
            ,"channel_32.peak_pos - channel_33.peak_pos > -30 && channel_32.peak_pos - channel_33.peak_pos < 30"
            "&& channel_32.amp < 60000 && channel_32.amp > 400"
            "&& channel_33.amp > 400 && channel_33.amp < 60000"
    #if UseNotEntangledPhotons
            "&& channel_32.peak_pos - channel_34.peak_pos > -60 && channel_32.peak_pos - channel_34.peak_pos < 60"
            "&& channel_34.amp > 400 && channel_34.amp < 60000"            
    #endif

            );
            true_sc_sum_hist_left->GetYaxis()->SetRangeUser(0,(Int_t)1.05*true_sc_sum_hist_left->GetMaximum());
            canv_00->cd(2);
            PMT_tree->Draw("channel_33.integral_in_gate >> TRUE_peak_scatterer_right_sum_spectrum"
            ,"channel_32.peak_pos - channel_33.peak_pos > -30 && channel_32.peak_pos - channel_33.peak_pos < 30"
            "&& channel_32.amp < 60000 && channel_32.amp > 400"
            "&& channel_33.amp > 400 && channel_33.amp < 60000"
    #if UseNotEntangledPhotons
            "&& channel_32.peak_pos - channel_34.peak_pos > -60 && channel_32.peak_pos - channel_34.peak_pos < 60"
            "&& channel_34.amp > 400 && channel_34.amp < 60000"            
    #endif
            );
            true_sc_sum_hist_right->GetYaxis()->SetRangeUser(0,(Int_t)1.05*true_sc_sum_hist_right->GetMaximum());

            canv_00->SaveAs(result_path + ".pdf(",".pdf");
            true_sc_sum_hist_left->Write();
            true_sc_sum_hist_right->Write();
*/
//////if hits NaI near photopeak
        for (Int_t channel_number = 0; channel_number < 32; channel_number++)
        {
            TString hernya = Form("integral_in_gate_scatterer_ch_cut_%i",channel_number);
            peak_histo[channel_number] = new TH1F (hernya.Data(),hernya.Data(),100,1,500);
        }


        for (Int_t iEvent = 0; iEvent < Events_for_cuts; iEvent++)
        {
            PMT_tree->GetEntry(iEvent);
            for (Int_t channel_number = 0; channel_number < 32; channel_number++)
            {
                Int_t sc_number = 0;
                if (channel_number < 16) sc_number = 32;
                if (channel_number >= 16) sc_number = 33;
                if (short_channel_info[channel_number].integral_in_gate > low_cut[channel_number]
                && short_channel_info[channel_number].integral_in_gate < high_cut[channel_number]

                && (short_channel_info[channel_number].peak_pos - short_channel_info[sc_number].peak_pos) > 100
                && (short_channel_info[channel_number].peak_pos - short_channel_info[sc_number].peak_pos) < 200

                && short_channel_info[channel_number].amp > 400 && short_channel_info[channel_number].amp < 60000
                && short_channel_info[sc_number].amp > 400 && short_channel_info[sc_number].amp < 60000                
         

                &&   short_channel_info[channel_number].peak_pos > 1600
                &&   short_channel_info[33].peak_pos > 1600
                &&   short_channel_info[32].peak_pos > 1600

                && short_channel_info[33].amp < 60000
                && short_channel_info[33].amp > 400
                && short_channel_info[32].amp < 60000     
                && short_channel_info[32].amp > 400 
                && short_channel_info[channel_number].amp > 400 && short_channel_info[channel_number].amp < 60000


#if UseNotEntangledPhotons

                && short_channel_info[34].amp < 50000
                && short_channel_info[34].amp > 400
                && (short_channel_info[34].peak_pos - short_channel_info[32].peak_pos) < 60
                && (short_channel_info[34].peak_pos - short_channel_info[32].peak_pos) >-30                 
#else 
                && short_channel_info[34].amp < 400
                && abs(short_channel_info[34].peak_pos - short_channel_info[32].peak_pos) > 100
#endif
      



                ) 
                    peak_histo[channel_number]->Fill(short_channel_info[sc_number].integral_in_gate);
            }
        }

    TString left_hist_name = "integral_in_gate_ch_32";
    TString right_hist_name = "integral_in_gate_ch_33";
    TH1F *right_hist = new TH1F (right_hist_name.Data(),right_hist_name.Data(),100,1,500);
    TH1F *left_hist = new TH1F (left_hist_name.Data(),left_hist_name.Data(),100,1,500);


            for (Int_t channel_number = 0; channel_number < 32; channel_number++)
            {
                cout << "channel_number = "<<channel_number <<endl;

				TF1 *gaus_fit = new TF1 ("gaus_func","gaus",150,350);
				peak_histo[channel_number]->Fit("gaus_func","R");

				TF1 *gaus_2_fit = new TF1 ("gaus_func_2","gaus",
				(gaus_fit->GetParameter(1))-gaus_fit->GetParameter(2),
				(gaus_fit->GetParameter(1))+gaus_fit->GetParameter(2));
				peak_histo[channel_number]->Fit("gaus_func_2","R");

                average_scatterer_peak_position[channel_number] = gaus_2_fit->GetParameter(1);
                sigma_scatterer[channel_number] = gaus_2_fit->GetParameter(2);

            if (channel_number < 16) left_hist->Add(peak_histo[channel_number]);
            if (channel_number >= 16) right_hist->Add(peak_histo[channel_number]);

/*                TCanvas *canv = new TCanvas("canv","canv");
                canv->cd();
                peak_histo[channel_number]->Draw();

                canv->SaveAs(source_path+"result_ev_b_ev.pdf(",".pdf");

                delete canv;
*/
                delete peak_histo[channel_number];

            }
                hist_like_ivashkin_wants_it(left_hist,"Energy [keV]", "Counts");
                hist_like_ivashkin_wants_it(right_hist,"Energy [keV]", "Counts");

                left_hist->Write();
                right_hist->Write();

                #if DrawToPDF
                TCanvas *canv = new TCanvas("canv","canv");
                canv->Divide(2);
                canv-> cd(1);
                left_hist->Draw();
                canv->cd(2);
                right_hist->Draw();

                canv->SaveAs(result_path+".pdf(",".pdf");
                delete canv;

                #endif

                delete right_hist;
                delete left_hist;
                


//////////////////////////////////sum_energy_histos
///////////////////////////////////

            for (Int_t channel_number = 0; channel_number < 32; channel_number++)
            {
                sum_energy[channel_number] = 510;
                
                sigma_energy[channel_number] = sqrt(0.5*(pow((float)sigma_scatterer[channel_number],2)+pow((float)sigma_i[channel_number],2)));
                cout<< channel_number<<endl;
                cout<< sigma_scatterer[channel_number] << " = sigma_scatterer" <<endl; 
                cout<< sigma_energy[channel_number] << " = sigma_energy" <<endl; 

                low_energy_cut[channel_number] = average_scatterer_peak_position[channel_number] - interval_width*sigma_scatterer[channel_number];
                low_energy_sum[channel_number] = sum_energy[channel_number] - interval_width*sigma_energy[channel_number];
                
                cout<< low_energy_cut[channel_number] << " = low_energy_cut" <<endl; 
                cout<< low_energy_sum[channel_number] << " = low_energy_sum" <<endl<<endl; 

                high_energy_cut[channel_number] = average_scatterer_peak_position[channel_number] + interval_width*sigma_scatterer[channel_number];
                high_energy_sum[channel_number] = sum_energy[channel_number] + interval_width*sigma_energy[channel_number];
                
                cout<< high_energy_cut[channel_number] << " = high_energy_cut" <<endl; 
                cout<< high_energy_sum[channel_number] << " = high_energy_sum" <<endl<<endl;           
            }
////////////////////////////////////
//////////////////////////////////////////////////
#if DrawTime
        for (Int_t channel_number = 0; channel_number < 35; channel_number++)
        {
                Int_t sc_number = 33;
                if (channel_number == 33) sc_number = 34;
                if (channel_number < 16 || channel_number == 34) sc_number = 32;
            TString hernya = Form("Time_difference_scatterer_%i_and_channel_%i_with_energy_cut",sc_number,channel_number);

            peak_histo[channel_number] = new TH1F (hernya.Data(),hernya.Data(),601,-300,300);

            if (channel_number < 32)
            {


                TString fignya = hernya + "without_energy_cuts";
                peak_histo_w_o_cuts[channel_number] = new TH1F (fignya.Data(),fignya.Data(),601,-300,300);                
            }

        }
        for (Int_t iEvent = 0; iEvent < Events_for_cuts; iEvent++)
        {
            PMT_tree->GetEntry(iEvent);
            for (Int_t channel_number = 0; channel_number < 35; channel_number++)
            {
                Int_t sc_number = 33;

                if (channel_number < 16 || channel_number == 34) sc_number = 32;

                if (channel_number < 32)
                {
                    if (short_channel_info[33].amp < 60000
                    && short_channel_info[33].amp > 400
                    && short_channel_info[32].amp < 60000     
                    && short_channel_info[32].amp > 400

                &&   short_channel_info[channel_number].peak_pos > 1600
                &&   short_channel_info[33].peak_pos > 1600
                &&   short_channel_info[32].peak_pos > 1600


#if UseNotEntangledPhotons

                && short_channel_info[34].amp < 50000
                && short_channel_info[34].amp > 400
                && (short_channel_info[34].peak_pos - short_channel_info[32].peak_pos) < 60
                && (short_channel_info[34].peak_pos - short_channel_info[32].peak_pos) >-30 

#else 
                && short_channel_info[34].amp < 400
                && abs(short_channel_info[34].peak_pos - short_channel_info[32].peak_pos) > 100
#endif

                    && short_channel_info[channel_number].amp > 400 && short_channel_info[channel_number].amp < 60000
                    && short_channel_info[sc_number].amp > 400 && short_channel_info[sc_number].amp < 60000                
                    )
                    {
                        peak_histo_w_o_cuts[channel_number]->Fill(short_channel_info[sc_number].peak_pos
                        -short_channel_info[channel_number].peak_pos);
                        if (short_channel_info[channel_number].integral_in_gate > low_cut[channel_number]
                        && short_channel_info[channel_number].integral_in_gate < high_cut[channel_number] )
                            peak_histo[channel_number]->Fill(short_channel_info[sc_number].peak_pos - short_channel_info[channel_number].peak_pos);
                    }                    
                }



                else
                {
                    if (channel_number == 33) sc_number = 34;
                if (
                     short_channel_info[33].amp < 60000
                    && short_channel_info[33].amp > 400
                    && short_channel_info[32].amp < 60000     
                    && short_channel_info[32].amp > 400


#if UseNotEntangledPhotons

                && short_channel_info[34].amp < 50000
                && short_channel_info[34].amp > 400
#else 
                && short_channel_info[34].amp < 400
                && abs(short_channel_info[34].peak_pos - short_channel_info[32].peak_pos) > 100
#endif


                    && short_channel_info[channel_number].amp > 400 && short_channel_info[channel_number].amp < 60000
                    )
                    {
                        if (channel_number!=34  
#if UseNotEntangledPhotons                        
                            && (short_channel_info[34].peak_pos - short_channel_info[32].peak_pos) < 60
                            && (short_channel_info[34].peak_pos - short_channel_info[32].peak_pos) >-30
#endif                            
                             )

                            peak_histo[channel_number]->Fill(short_channel_info[sc_number].peak_pos
                            -short_channel_info[channel_number].peak_pos);
                        if (channel_number==34)
                            peak_histo[channel_number]->Fill(short_channel_info[sc_number].peak_pos
                            -short_channel_info[channel_number].peak_pos);                        
                    }
                }


            }
        }


            for (Int_t channel_number = 0; channel_number < 35; channel_number++)
            {
                cout << "channel_number = "<<channel_number <<endl;
#if DrawToPDF
                TCanvas *canv = new TCanvas("canv","canv");
#endif
                if (channel_number < 32)
                {
                    hist_like_ivashkin_wants_it(peak_histo_w_o_cuts[channel_number],"Time [ns]", "Counts");
                    hist_like_ivashkin_wants_it(peak_histo[channel_number],"Time [ns]", "Counts");

                    peak_histo_w_o_cuts[channel_number]->Write();
                    peak_histo[channel_number]->Write();

                    #if DrawToPDF
                    canv->Divide(2);
                    canv->cd(1);
                    peak_histo_w_o_cuts[channel_number]->Draw();
                    canv->cd(2);
                    peak_histo[channel_number]->Draw();
                    #endif
                }
                else
                {
                    hist_like_ivashkin_wants_it(peak_histo[channel_number],"Time [ns]", "Counts");

                    #if DrawToPDF
                    canv->cd();
                    peak_histo[channel_number]->Draw();
                    #endif
                    peak_histo[channel_number]->Write();

                }
                #if DrawToPDF
                canv->SaveAs(result_path+".pdf(",".pdf");
                delete canv;
                #endif

                delete peak_histo[channel_number];
            }
#endif
///////////////////////////////////
        for (Int_t channel_number = 0; channel_number < 32; channel_number++)
        {
            TString hernya = Form("integral_in_gate_sum_with_scatterer_ch_cut_%i",channel_number);
//            peak_histo[channel_number] = new TH1F (hernya.Data(),hernya.Data(),70,-350,350);
            peak_histo[channel_number] = new TH1F (hernya.Data(),hernya.Data(),70,200,900);

        }

        for (Int_t iEvent = 0; iEvent < Events_for_cuts; iEvent++)
        {
            PMT_tree->GetEntry(iEvent);
            for (Int_t channel_number = 0; channel_number < 32; channel_number++)
            {
                Int_t sc_number = 0; Short_t time_min = 0; Short_t time_max = 0;
                if (channel_number < 16) sc_number = 32;
                if (channel_number >= 16) sc_number = 33;
                if (
                    short_channel_info[channel_number].integral_in_gate > low_cut[channel_number]
                && short_channel_info[channel_number].integral_in_gate < high_cut[channel_number]

                && (short_channel_info[channel_number].peak_pos - short_channel_info[sc_number].peak_pos) > 100
                && (short_channel_info[channel_number].peak_pos - short_channel_info[sc_number].peak_pos) < 200
                && abs(short_channel_info[33].peak_pos - short_channel_info[32].peak_pos) < 30

                &&   short_channel_info[channel_number].peak_pos > 1600
                &&   short_channel_info[33].peak_pos > 1600
                &&   short_channel_info[32].peak_pos > 1600

#if UseNotEntangledPhotons

                && short_channel_info[34].amp < 50000
                && short_channel_info[34].amp > 400
                && (short_channel_info[34].peak_pos - short_channel_info[32].peak_pos) < 60
                && (short_channel_info[34].peak_pos - short_channel_info[32].peak_pos) >-30 

               // && short_channel_info[34].integral_in_gate >= 15
                && short_channel_info[34].integral_in_gate < 150
#else 
                && short_channel_info[34].amp < 400
                && abs(short_channel_info[34].peak_pos - short_channel_info[32].peak_pos) > 100
#endif
                /*&& short_channel_info[33].amp < 60000
                && short_channel_info[33].amp > 400
                && short_channel_info[32].amp < 60000     
                && short_channel_info[32].amp > 400*/

                && short_channel_info[sc_number].integral_in_gate > low_energy_cut[channel_number]
                && short_channel_info[sc_number].integral_in_gate < high_energy_cut[channel_number]

                && short_channel_info[channel_number].amp > 400 && short_channel_info[channel_number].amp < 60000
                && short_channel_info[sc_number].amp > 400 && short_channel_info[sc_number].amp < 60000                
                )
                {
                    peak_histo[channel_number]->Fill(
                        //510-
                        (short_channel_info[sc_number].integral_in_gate+short_channel_info[channel_number].integral_in_gate)
                    //+abs((Float_t)(33-sc_number)*short_channel_info[34].integral_in_gate)
                    )
                    ;        

                }

            }
        }


            for (Int_t channel_number = 0; channel_number < 32; channel_number++)
            {
                cout << "channel_number = "<<channel_number <<endl;

                hist_like_ivashkin_wants_it(peak_histo[channel_number],"Energy [keV]", "Counts");
                peak_histo[channel_number]->Write();

                #if DrawToPDF
                TCanvas *canv = new TCanvas("canv","canv");

                canv->cd();
                peak_histo[channel_number]->Draw();
                canv->SaveAs(result_path+".pdf(",".pdf");
                delete canv;

                #endif

#if TotalEnergyCut
				TF1 *gaus_fit = new TF1 ("gaus_func","gaus",250,700);
				peak_histo[channel_number]->Fit("gaus_func","R");

				TF1 *gaus_2_fit = new TF1 ("gaus_func_2","gaus",
				(gaus_fit->GetParameter(1))-gaus_fit->GetParameter(2),
				(gaus_fit->GetParameter(1))+gaus_fit->GetParameter(2));
				peak_histo[channel_number]->Fit("gaus_func_2","R");


				total_low_cut[channel_number] = gaus_2_fit->GetParameter(1)-0.38*gaus_2_fit->GetParameter(2);
				total_high_cut[channel_number] = gaus_2_fit->GetParameter(1)+0.38*gaus_2_fit->GetParameter(2);
#endif 



                delete peak_histo[channel_number];
            }


///////////////////////////////////
///////////////////////////////////
#if CalculateRatio
    for (Int_t NumEvent = 0; NumEvent < calculate_Events_Number; NumEvent++)
    {
        PMT_tree->GetEntry(NumEvent);
        if (NumEvent%1000000 == 0) cout << (float)NumEvent/calculate_Events_Number*100<<"%"<<endl;
        Short_t counter_1 = 0;
        Short_t counter_2 = 0;
        
        for (Int_t channel_number = 0; channel_number < 16; channel_number++) if (short_channel_info[channel_number].amp > 400) counter_1++;
        for (Int_t channel_number = 16; channel_number < 32; channel_number++) if (short_channel_info[channel_number].amp > 400) counter_2++;


        for (Int_t channel_number = 0; channel_number < 16; channel_number++)
        {
            for (Int_t channel_number_2 = 16; channel_number_2 < 32; channel_number_2++)
            {
                if(
                    abs(short_channel_info[33].peak_pos - short_channel_info[32].peak_pos) < 30

                    && counter_1 ==1 && counter_2 ==1
           
#if UseIntegralCut
                && short_channel_info[channel_number].integral_in_gate > low_cut[channel_number]
                && short_channel_info[channel_number].integral_in_gate < high_cut[channel_number]
                && short_channel_info[channel_number_2].integral_in_gate > low_cut[channel_number_2]
                && short_channel_info[channel_number_2].integral_in_gate < high_cut[channel_number_2]

                && short_channel_info[32].integral_in_gate > low_energy_cut[channel_number]
                && short_channel_info[32].integral_in_gate < high_energy_cut[channel_number]
                && short_channel_info[33].integral_in_gate > low_energy_cut[channel_number_2]
                && short_channel_info[33].integral_in_gate < high_energy_cut[channel_number_2]                
#endif                
                && abs(short_channel_info[33].peak_pos - short_channel_info[32].peak_pos) < 30


                &&   short_channel_info[channel_number].peak_pos > 1600
                &&   short_channel_info[33].peak_pos > 1600
                &&   short_channel_info[32].peak_pos > 1600

                && short_channel_info[33].amp < 60000
                && short_channel_info[33].amp > 400
                && short_channel_info[32].amp < 60000     
                && short_channel_info[32].amp > 400 
                && short_channel_info[channel_number].amp > 400 && short_channel_info[channel_number].amp < 60000
                && short_channel_info[channel_number_2].amp > 400 && short_channel_info[channel_number_2].amp < 60000

                && (short_channel_info[channel_number].peak_pos - short_channel_info[32].peak_pos) > 100
                && (short_channel_info[channel_number].peak_pos - short_channel_info[32].peak_pos) < 200
                && (short_channel_info[channel_number_2].peak_pos - short_channel_info[33].peak_pos) > 100
                && (short_channel_info[channel_number_2].peak_pos - short_channel_info[33].peak_pos) < 200

#if UseNotEntangledPhotons

                && short_channel_info[34].amp < 50000
                && short_channel_info[34].amp > 400
                && (short_channel_info[34].peak_pos - short_channel_info[32].peak_pos) < 60
                && (short_channel_info[34].peak_pos - short_channel_info[32].peak_pos) >-30                 
                //&& short_channel_info[34].integral_in_gate >= 15
                && short_channel_info[34].integral_in_gate < 150

#else 
                && short_channel_info[34].amp < 400
                && abs(short_channel_info[34].peak_pos - short_channel_info[32].peak_pos) > 100
#endif

#if TotalEnergyCut
                //&& short_channel_info[32].integral_in_gate + short_channel_info[channel_number].integral_in_gate  >  total_high_cut[channel_number]
                && short_channel_info[32].integral_in_gate + short_channel_info[channel_number].integral_in_gate < total_low_cut[channel_number]

#endif
                )
                {
                    hist_charge[channel_number][channel_number_2-16]->Fill(
                    short_channel_info[32].integral_in_gate + short_channel_info[channel_number].integral_in_gate);
                }
            }
        }
    }

        for (Int_t channel_number = 0; channel_number < 16; channel_number++)
        {
            for(Int_t channel_number_2 = 16; channel_number_2 < 32; channel_number_2++)
            {
                NumEvents[channel_number][channel_number_2-16] = hist_charge[channel_number][channel_number_2-16]->GetEntries();
                /*TCanvas *canv = new TCanvas("canv","canv");
                canv->cd();
                hist_charge[channel_number][channel_number_2-16]->Draw();
                canv->SaveAs(source_path+"result_ev_b_ev.pdf(",".pdf");*/
                Int_t angle = channel_number_2-channel_number-16;
                if (angle < 0) angle += 16;
                total_events_for_angle_diff[angle]+=NumEvents[channel_number][channel_number_2-16];

            }
        }

/////////////////////
/////////////////////CHSH inequality
            #if calculate_CHSH
                const Int_t different_angles = 8;
                const Int_t total_angles = different_angles*2;
                Float_t phi_angle[total_angles] = {0.};
                for (Int_t i = 0; i < different_angles; i++) phi_angle[i] = (float)(i+1)*22.5;
                for (Int_t i = 0; i < different_angles; i++) phi_angle[i+different_angles] = -(float)(i+1)*22.5;



                Float_t CHSH[total_angles] = {0};
                Float_t CHSH_error[total_angles] = {0};
                Float_t global_CHSH[total_angles] = {0};
                Float_t global_CHSH_error[total_angles] = {0};
                for (Int_t a_angle = 0; a_angle < 16; a_angle++)
                {
                    Float_t CHSH_error_local[total_angles] = {0};
                    Float_t CHSH_local[total_angles] = {0};
                    for (Int_t angle_between_counters = 1; angle_between_counters <= different_angles; angle_between_counters++)
                    {
                        
                        CHSH_local[angle_between_counters-1] = 
                        calculate_local_CHSH(NumEvents, a_angle, angle_between_counters, "clockwise");
                        CHSH_local[angle_between_counters+different_angles-1] = calculate_local_CHSH(NumEvents, a_angle, angle_between_counters, "counterclockwise");


                        CHSH_error_local[angle_between_counters-1] = 
                        calculate_local_CHSH_error(NumEvents, a_angle, angle_between_counters, "clockwise");
                        CHSH_error_local[angle_between_counters-1+different_angles] = calculate_local_CHSH_error(NumEvents, a_angle, angle_between_counters, "counterclockwise");


                        CHSH[angle_between_counters-1] += CHSH_local[angle_between_counters-1];
                        CHSH_error[angle_between_counters-1] += CHSH_error_local[angle_between_counters-1];       
                        CHSH_error_local[angle_between_counters-1] = sqrt(CHSH_error_local[angle_between_counters-1]);           
                    }

                        TCanvas *canvas_CHSH = new TCanvas ( Form("canvas_CHSH_a=%i",a_angle), Form("canvas_CHSH_a=%i",a_angle));
                        canvas_CHSH->cd();
                        TGraphErrors * gr_CHSH = new TGraphErrors(total_angles,phi_angle,CHSH_local, angle_arr_err, CHSH_error_local);
                        TF1 *my_fit = new TF1 ("my_fit","[0]*(cos(6*3.141593/180*x)-3*cos(2*3.141593/180*x))",-190,190);
                        gr_CHSH->Fit("my_fit","R");
                        Float_t a_par_CHSH = abs(my_fit->GetParameter(0));
                        Float_t a_error_CHSH = abs(my_fit->GetParError(0));
                        graph_like_ivashkin_wants_it(gr_CHSH,"angle [degrees]","S", Form("canvas_CHSH_a=%i <> E(90) = %4.3f+-%4.3f",a_angle,E_coeff(NumEvents,a_angle,4, "clockwise"),sqrt(sqr_E_error(NumEvents,a_angle,4, "clockwise"))),1);
                        gr_CHSH->Draw("AP");
                        canvas_CHSH->SaveAs(result_path+".pdf",".pdf");

                }

                cout <<endl<<endl<<endl<<endl;
                for (Int_t i = 0; i < total_angles; i++) 
                {
                    CHSH[i] = CHSH[i]/16;
                    CHSH_error[i] = sqrt(CHSH_error[i])/16;
                }
                for (Int_t i = 0; i < total_angles; i++) cout<< CHSH[i] <<" ";
                cout <<endl<<endl<<endl<<endl; 
                for (Int_t i = 0; i < total_angles; i++) cout<< CHSH_error[i] <<" ";
                cout <<endl<<endl<<endl<<endl;                 

                cout <<endl<<endl<<endl<<endl;
                for (Int_t angle_between_counters = 1; angle_between_counters <= different_angles; angle_between_counters++)
                {
                    global_CHSH[angle_between_counters-1] = global_calculate_CHSH(angle_between_counters, NumEvents, "clockwise");
                    global_CHSH_error[angle_between_counters-1] = sqrt(global_calculate_CHSH_error(angle_between_counters, NumEvents, "clockwise"));
                    cout<< global_CHSH_error[angle_between_counters-1] <<" ";

                }
                cout <<endl<<endl<<endl<<endl;
            #endif

///////////////////////////////////////
///////////////////////////////////////

#endif    

#if calculate_CHSH
Float_t average_E = 0; Float_t sigma_E_average = 0;

    for (Int_t a = 0; a < 16; a++)
    {

        average_E += E_coeff(NumEvents,a,4, "clockwise")/16.;
        sigma_E_average += sqr_E_error(NumEvents,a,4, "clockwise");
    }
    sigma_E_average = sqrt(sigma_E_average)/16;

    //////////////////Draw average CHSH
    TGraphErrors * gr_CHSH = new TGraphErrors(different_angles,phi_angle,CHSH, angle_arr_err, CHSH_error);
	TF1 *my_fit = new TF1 ("my_fit","[0]*(cos(6*3.141593/180*x)-3*cos(2*3.141593/180*x))",10,190);
	gr_CHSH->Fit("my_fit","R");
    Float_t a_par_CHSH = abs(my_fit->GetParameter(0));
    Float_t a_error_CHSH = abs(my_fit->GetParError(0));
    TCanvas *canvas_CHSH = new TCanvas ( "canvas_CHSH", "canvas_CHSH");
    canvas_CHSH->cd();
    graph_like_ivashkin_wants_it(gr_CHSH,"angle [degrees]","S", Form("Average E <> E(90) = %4.3f+-%4.3f",average_E,sigma_E_average),1);
    gr_CHSH->Draw("AP");
    canvas_CHSH->SaveAs(result_path+".pdf",".pdf");


    ///////////////////Draw global CHSH
	TF1 *global_my_fit = new TF1 ("global_my_fit","[0]*(cos(6*3.141593/180*x)-3*cos(2*3.141593/180*x))",10,190);
    TGraphErrors * global_gr_CHSH = new TGraphErrors(different_angles,phi_angle,global_CHSH, angle_arr_err, global_CHSH_error);
	global_gr_CHSH->Fit("global_my_fit","R");
    a_par_CHSH = abs(global_my_fit->GetParameter(0));
    a_error_CHSH = abs(global_my_fit->GetParError(0));
    TCanvas *global_canvas_CHSH = new TCanvas ( "global_canvas_CHSH", "global_canvas_CHSH");
    global_canvas_CHSH->cd();
    graph_like_ivashkin_wants_it(global_gr_CHSH,"angle [degrees]","S", Form("Sum N <> E(90) = %4.3f+-%4.3f",global_E_coeff(4,NumEvents, "clockwise"), sqrt(global_sqr_E_error(4,NumEvents, "clockwise"))),1);
    global_gr_CHSH->Draw("AP");
    global_canvas_CHSH->SaveAs(result_path+".pdf",".pdf");


   

#endif

    for (Int_t i = 0; i < 16; i++) total_events_for_angle_diff_err[i] = pow(total_events_for_angle_diff[i],0.5);

    TGraphErrors * gr = new TGraphErrors(16,angle_arr,total_events_for_angle_diff, angle_arr_err, total_events_for_angle_diff_err);
	TF1 *sin_fit = new TF1 ("sin_fit","[0]+[1]*cos(3.14159265359/180*2*x)",0,360);
	gr->Fit("sin_fit","R");
    Float_t a_par = abs(sin_fit->GetParameter(0));
    Float_t b_par = abs(sin_fit->GetParameter(1));
    Float_t a_error = abs(sin_fit->GetParError(0));
    Float_t b_error = abs(sin_fit->GetParError(1));
    Float_t diff = (a_par+b_par)/(a_par-b_par);
    //Float_t diff_err = diff*sqrt(2*pow((a_error/a_par),2)+2*pow((b_error/b_par),2));
    Float_t diff_err = 2*a_par/(pow(a_par-b_par,2))*sqrt(pow((a_error*b_par/a_par),2)+pow((b_error),2));
    TCanvas *canvas = new TCanvas ( "canvas", "canvas");
    canvas->cd();
    graph_like_ivashkin_wants_it(gr,"angle [degrees]","Counts", Form("%5.1f-%5.1f*cos(2x) ratio = %4.3f +- %4.3f",a_par,b_par,diff,diff_err),1);
    gr->Draw("AP");

    canvas->SaveAs(result_path+".pdf)",".pdf");


    result_root->Close();

}



//////////////////////
