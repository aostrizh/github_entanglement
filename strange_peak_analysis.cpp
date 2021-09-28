#include <iostream>
#include <fstream>
#include <TH1F.h>
#include <TH2F.h>

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

#define UseManyRoots 0
#define PresentDiffuser 1
#define CalculateRatio 1
#define DrawTime 1
#define UseIntegralCut 1
#define UseNotEntangledPhotons 1
#define TotalEnergyCut 0
#define DrawFalsePeak 0
void strange_peak_analysis()
{
	gStyle->SetOptFit(1);
	//gStyle->SetOptStat(1111);
	TString source_path = "/home/doc/entanglement/root_files_data/with_scatterer/big_file/";

    TString result_path = source_path + "result_with_2d_scatterer_arm";
    TFile *result_root = new TFile (result_path+".root", "RECREATE");

////////////////
    Int_t low_cut[32] = {0}; Int_t high_cut[32] = {0}; Float_t total_low_cut[32] = {0}; Float_t total_high_cut[32] = {0};
    Float_t sigma_i[32] = {0}; Float_t mean_int[32] = {0};

    Int_t low_amp_cut[32] = {0};
    Int_t high_amp_cut[32] = {0};

    Float_t interval_width = 1.4;
    Float_t average_scatterer_peak_position[32] = {0};
    Float_t sigma_scatterer[32] = {0};
    Float_t low_energy_cut[32] = {0};
    Float_t high_energy_cut[32] = {0};

    Float_t sum_energy[32] = {0};
    Float_t sigma_energy[32] = {0};
    Float_t low_energy_sum[32] = {0};
    Float_t high_energy_sum[32] = {0};


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
    const Int_t calculate_Events_Number =  
    //50000000
    PMT_tree->GetEntries()/1
    ;
    const Int_t events_wanted = 10000000;
    Int_t Events_for_cuts = 0;


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
        const Int_t histos_bins_number = 100;
        const Int_t histos_left_boarder = 5;
        const Int_t histos_right_boarder = 450;
        TH1F *peak_histo_w_o_cuts[35];

#if DrawFalsePeak       
        TH1F *peak_histo[35]; 
        TH1F *scatterer_histo[35];
        TH1F *difuser_histo[35];
        for (Int_t channel_number = 0; channel_number < 32; channel_number++)
        {
                Int_t sc_number = 33;
                if (channel_number < 16) sc_number = 32;
            TString hernya = Form("integral_in_gate_channel_%i",channel_number);
            peak_histo[channel_number] = new TH1F (hernya.Data(),hernya.Data(),histos_bins_number,histos_left_boarder,histos_right_boarder);
            TString hernya1 = Form("integral_in_gate_in_sc_channel_cut_%i",channel_number);
            scatterer_histo[channel_number] = new TH1F (hernya1.Data(),hernya1.Data(),histos_bins_number,histos_left_boarder,histos_right_boarder);
            TString hernya2 = Form("integral_in_gate_in_difuser_channel_cut_%i",channel_number);
            difuser_histo[channel_number] = new TH1F (hernya2.Data(),hernya2.Data(),histos_bins_number,histos_left_boarder,histos_right_boarder);
        
        }
#endif


        TH1F *peak_histo_true[35];
        TH1F *scatterer_histo_true[35];
        TH1F *difuser_histo_true[35];
        TH2F *peak2d_true[35]; 
        TH2F *sc_vs_diff[35];
        TH2F *left_arm_vs_diff[35];


        for (Int_t channel_number = 0; channel_number < 32; channel_number++)
        {
                Int_t sc_number = 33;
                if (channel_number < 16) sc_number = 32;
            TString hernya = Form("TRUE_peak_integral_in_gate_channel_%i",channel_number);
            peak_histo_true[channel_number] = new TH1F (hernya.Data(),hernya.Data(),histos_bins_number,histos_left_boarder,histos_right_boarder);
            TString hernya1 = Form("TRUE_peak_integral_in_gate_in_sc_channel_cut_%i",channel_number);
            scatterer_histo_true[channel_number] = new TH1F (hernya1.Data(),hernya1.Data(),histos_bins_number,histos_left_boarder,histos_right_boarder);
            TString hernya2 = Form("TRUE_peak_integral_in_gate_in_difuser_channel_cut_%i",channel_number);
            difuser_histo_true[channel_number] = new TH1F (hernya2.Data(),hernya2.Data(),histos_bins_number,histos_left_boarder,histos_right_boarder);
            TString hernya3 = Form("TRUE_peak_integral_in_gate_ch_%i_vs_sc",channel_number);
            peak2d_true[channel_number] = new TH2F (hernya3.Data(),hernya3.Data(),100,0,500,100,0,500);
            TString hernya4 = Form("left_scatterer_vs_weak_scatterer_ch_%i",channel_number);
            sc_vs_diff[channel_number] = new TH2F (hernya4.Data(),hernya4.Data(),100,0,500,100,0,300);
            TString hernya5 = Form("left_arm_vs_weak_scatterer_ch_%i",channel_number);
            left_arm_vs_diff[channel_number] = new TH2F (hernya5.Data(),hernya5.Data(),100,0,500,100,0,300);
        
        }

        for (Int_t iEvent = 0; iEvent < calculate_Events_Number; iEvent++)
        {
            PMT_tree->GetEntry(iEvent);
            for (Int_t channel_number = 0; channel_number < 32; channel_number++)
            {
                Int_t sc_number = 33;
                if (channel_number < 16) sc_number = 32;
                if (

//////////////////

                abs(short_channel_info[33].peak_pos - short_channel_info[32].peak_pos) < 30

                &&   short_channel_info[channel_number].peak_pos > 1600
                &&   short_channel_info[33].peak_pos > 1600
                &&   short_channel_info[32].peak_pos > 1600

                && short_channel_info[33].amp < 60000
                && short_channel_info[33].amp > 400
                && short_channel_info[32].amp < 60000     
                && short_channel_info[32].amp > 400 
                && short_channel_info[channel_number].amp > 400 && short_channel_info[channel_number].amp < 60000
                && short_channel_info[channel_number].integral_in_gate > 4


#if UseIntegralCut

                && short_channel_info[32].integral_in_gate > 10
                && short_channel_info[32].integral_in_gate < 500
                && short_channel_info[33].integral_in_gate > 10
                && short_channel_info[33].integral_in_gate < 500

#endif         
#if UseNotEntangledPhotons

                && short_channel_info[34].amp < 50000
                && short_channel_info[34].amp > 400
                && (short_channel_info[34].peak_pos - short_channel_info[32].peak_pos) < 60
                && (short_channel_info[34].peak_pos - short_channel_info[32].peak_pos) >-60                 
#else 
                && short_channel_info[34].amp < 400
                && abs(short_channel_info[34].peak_pos - short_channel_info[32].peak_pos) > 100
#endif

                ) 
                    {
                        //cout<< "TRUE"<< endl;
#if DrawFalsePeak
                        if (abs(short_channel_info[channel_number].peak_pos - short_channel_info[33].peak_pos) < 80)
                        {
                            peak_histo[channel_number]->Fill(short_channel_info[channel_number].integral_in_gate);
                            scatterer_histo[channel_number]->Fill(short_channel_info[sc_number].integral_in_gate);
                            difuser_histo[channel_number]->Fill(short_channel_info[34].integral_in_gate);


                        }
#endif
                        if((short_channel_info[channel_number].peak_pos - short_channel_info[33].peak_pos) > 80
                            && (short_channel_info[channel_number].peak_pos - short_channel_info[33].peak_pos) < 200)
                        {
                            peak_histo_true[channel_number]->Fill(short_channel_info[channel_number].integral_in_gate);
                            scatterer_histo_true[channel_number]->Fill(short_channel_info[sc_number].integral_in_gate);
                            difuser_histo_true[channel_number]->Fill(short_channel_info[34].integral_in_gate);
                            peak2d_true[channel_number]->Fill(short_channel_info[sc_number].integral_in_gate,short_channel_info[channel_number].integral_in_gate);
                            left_arm_vs_diff[channel_number]->Fill(short_channel_info[channel_number].integral_in_gate,short_channel_info[34].integral_in_gate);
                            sc_vs_diff[channel_number]->Fill(short_channel_info[sc_number].integral_in_gate,short_channel_info[34].integral_in_gate);
                        }
                    }

            }
        }
#if DrawFalsePeak
        TH1F *sum_hist_left = new TH1F ("left_sum_spectrum","left_sum_spectrum",histos_bins_number,histos_left_boarder,histos_right_boarder);
        TH1F *sum_hist_right = new TH1F ("right_sum_spectrum","right_sum_spectrum",histos_bins_number,histos_left_boarder,histos_right_boarder);
        TH1F *sc_sum_hist_left = new TH1F ("scatterer_left_sum_spectrum","scatterer_left_sum_spectrum",histos_bins_number,histos_left_boarder,histos_right_boarder);
        TH1F *sc_sum_hist_right = new TH1F ("scatterer_right_sum_spectrum","scatterer_right_sum_spectrum",histos_bins_number,histos_left_boarder,histos_right_boarder);
        TH1F *difuser_sum_hist_left = new TH1F ("difuser_sum_spectrum","difuser_sum_spectrum",histos_bins_number,histos_left_boarder,histos_right_boarder);
#endif        
        TH1F *true_sum_hist_left = new TH1F ("TRUE_peak_left_sum_spectrum","TRUE_peak_left_sum_spectrum",histos_bins_number,histos_left_boarder,histos_right_boarder);
        TH1F *true_sum_hist_right = new TH1F ("TRUE_peak_right_sum_spectrum","TRUE_peak_right_sum_spectrum",histos_bins_number,histos_left_boarder,histos_right_boarder);
        TH1F *true_sc_sum_hist_left = new TH1F ("TRUE_peak_scatterer_left_sum_spectrum","TRUE_peak_scatterer_left_sum_spectrum",histos_bins_number,histos_left_boarder,histos_right_boarder);
        TH1F *true_sc_sum_hist_right = new TH1F ("TRUE_peak_scatterer_right_sum_spectrum","TRUE_peak_scatterer_right_sum_spectrum",histos_bins_number,histos_left_boarder,histos_right_boarder);
        TH1F *true_difuser_sum_hist_left = new TH1F ("TRUE_peak_difuser_sum_spectrum","TRUE_peak_difuser_sum_spectrum",histos_bins_number,histos_left_boarder,histos_right_boarder);
        TH2F *left_sumch_vs_scatterer = new TH2F ("TRUE_peak_integral_in_gate_left_vs_sc","TRUE_peak_integral_in_gate_left_vs_sc",100,0,500,100,0,500);         
        TH2F *right_sumch_vs_scatterer = new TH2F ("TRUE_peak_integral_in_gate_right_vs_sc","TRUE_peak_integral_in_gate_right_vs_sc",100,0,500,100,0,500);         
        TH2F *sc_vs_diffuser = new TH2F ("TRUE_peak_integral_in_gate_sc_weak_sc","TRUE_peak_integral_in_gate_sc_weak_sc",100,0,500,100,0,300);         
        TH2F *larm_vs_diffuser = new TH2F ("TRUE_peak_integral_in_gate_left_arm_weak_sc","TRUE_peak_integral_in_gate_left_arm_weak_sc",100,0,500,100,0,300);         
            
            
            for (Int_t channel_number = 0; channel_number < 32; channel_number++)
            {
				

                if ( channel_number < 16) 
                {
#if DrawFalsePeak
                   sum_hist_left->Add(peak_histo[channel_number]); 
                   sc_sum_hist_left->Add(scatterer_histo[channel_number]);
                   difuser_sum_hist_left->Add(difuser_histo[channel_number]);
#endif
                   true_sum_hist_left->Add(peak_histo_true[channel_number]);
                   true_sc_sum_hist_left->Add(scatterer_histo_true[channel_number]);
                   true_difuser_sum_hist_left->Add(difuser_histo_true[channel_number]);
                   left_sumch_vs_scatterer->Add(peak2d_true[channel_number]);
                   sc_vs_diffuser->Add(sc_vs_diff[channel_number]);
                   larm_vs_diffuser->Add(left_arm_vs_diff[channel_number]);
                }
                if ( channel_number >= 16) 
                {
#if DrawFalsePeak                    
                   sum_hist_right->Add(peak_histo[channel_number]); 
                   sc_sum_hist_right->Add(scatterer_histo[channel_number]);
#endif
                   true_sum_hist_right->Add(peak_histo_true[channel_number]);
                   true_sc_sum_hist_right->Add(scatterer_histo_true[channel_number]);
                   right_sumch_vs_scatterer->Add(peak2d_true[channel_number]);

                }
            }



    TCanvas *canvas = new TCanvas ( "canvas", "canvas");
    TCanvas *canvas_2 = new TCanvas ( "canvas_2", "canvas_2");
    TCanvas *canvas_3 = new TCanvas ( "canvas_3", "canvas_3");
    TCanvas *canvas_4 = new TCanvas ( "canvas_4", "canvas_4");
    

#if DrawFalsePeak
    canvas->Divide(2,2);
    canvas_2->Divide(2,2);
    canvas_3->Divide(2);
#else
    canvas->Divide(2);
    canvas_2->Divide(2);


#endif

    canvas->cd(1);
    true_sum_hist_left->Draw(); true_sum_hist_left->Write();
    canvas->cd(2);
    true_sum_hist_right->Draw();    true_sum_hist_right->Write();

    canvas_2->cd(1);
    true_sc_sum_hist_left->Draw();    true_sc_sum_hist_left->Write();

    canvas_2->cd(2);
    true_sc_sum_hist_right->Draw();    true_sc_sum_hist_right->Write();

#if DrawFalsePeak
    canvas_3->cd(1);
#else
    canvas_3 ->cd();
#endif
    true_difuser_sum_hist_left->Draw();    true_difuser_sum_hist_left->Write();

    canvas_4->Divide(2,2);
    canvas_4->cd(1);
    left_sumch_vs_scatterer->Draw("colz");    left_sumch_vs_scatterer->Write();

    canvas_4->cd(2);
    right_sumch_vs_scatterer->Draw("colz");    right_sumch_vs_scatterer->Write();

    canvas_4->cd(3);
    sc_vs_diffuser->Draw("colz");    sc_vs_diffuser->Write();

    canvas_4->cd(4);
    larm_vs_diffuser->Draw("colz"); larm_vs_diffuser->Write();;



#if DrawFalsePeak
    canvas->cd(3);
    sum_hist_left->Draw();
    canvas->cd(4);
    sum_hist_right->Draw();
    canvas_2->cd(3);
    sc_sum_hist_left->Draw();
    canvas_2->cd(4);
    sc_sum_hist_right->Draw();
    canvas_3->cd(2);
    difuser_sum_hist_left->Draw();
#endif

    canvas->SaveAs(result_path+".pdf(",".pdf");
    canvas_2->SaveAs(result_path+".pdf(",".pdf");    
    canvas_3->SaveAs(result_path+".pdf(",".pdf");
    canvas_4->SaveAs(result_path+".pdf)",".pdf");
    delete canvas;
    delete canvas_2;
    delete canvas_3;
    delete canvas_4;


}
