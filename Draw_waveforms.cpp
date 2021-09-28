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
#include <TGraph.h>
#include <TObjString.h>
#include "../calo_analysis/Readme_reader.h"
#include "ChannelEntry.h"
#define UsedScatterer 1

void Draw_waveforms()
{

    TString source_path = "/home/daq/entanglement/scatterer_4_8_cm_further_from_diffuser/";
    TChain *PMT_tree = new TChain;

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
			if(!file->IsDirectory() && fname.EndsWith("07a8de9a_20210312_235436.root")) 
			{
				files_names.Add(new TObjString(fname));
				PMT_tree->AddFile( (source_path + fname + "/adc64_data").Data() );
			}
		}
	}
    Int_t total_entries = PMT_tree->GetEntries();

#if UsedScatterer
    const Int_t total_channels = 35;
#else 
    const Int_t total_channels = 34;
#endif

     Int_t looK_channel = 5;
    const Int_t tot_canv = 10;
    //Short_t amp;
    //Short_t wf_size;
    //Short_t wf[MAX_N_SAMPLES];

    Int_t global_counter_1 = 0;
    Int_t global_counter_2 = 0;
    Int_t global_counter_3 = 0;
    Int_t samples[2050];

    ChannelEntry *channel_info = new ChannelEntry[total_channels];
    for(Int_t ch = 0; ch < total_channels; ch++)
	    (channel_info+ch)->SetBranch(PMT_tree, ch);

    for (Int_t i = 0; i < 2050; i++) samples[i] = i;
        for (Int_t entryNum = 0; entryNum < total_entries; entryNum++)
        {
            PMT_tree->GetEntry(entryNum);
            for (Int_t channel_number = 32; channel_number < 33; channel_number++)
            {
                looK_channel = channel_number;
            
            Int_t sc_channel = 33;
            if (looK_channel < 16) sc_channel = 32;
            const Int_t waveform_size = channel_info[looK_channel].wf_size;
              // cout << global_counter_1 <<" "<<global_counter_2<<" "<< global_counter_3<<" "
             //  <<abs(channel_info[sc_channel].peak_pos-channel_info[looK_channel].peak_pos)<<" "<<entryNum<<" "<<channel_info[looK_channel].wf_size<< endl;

            Short_t flag = 0;
            if ((global_counter_2 < tot_canv 
            //|| global_counter_1 < tot_canv || global_counter_3 < tot_canv
            ))
            {
            if (waveform_size == 350)
            {
    
                
                TCanvas *canv = new TCanvas("canv","canv");
               //cout << global_counter_1 <<" "<<global_counter_2<<" "<< global_counter_3<<" "                <<abs(channel_info[sc_channel].peak_pos-channel_info[looK_channel].peak_pos)<<" "<<entryNum<< endl;

                Int_t wf[waveform_size];
                cout << global_counter_1 <<" "<<global_counter_2<<" "<< global_counter_3<<" "
                <<abs(channel_info[sc_channel].peak_pos-channel_info[looK_channel].peak_pos)
                <<" "<<entryNum<<" "<<channel_info[looK_channel].wf_size<< endl;

                for (Int_t i = 0; i < waveform_size; i++)
                    wf[i] = channel_info[looK_channel].wf[i];

/*            if (global_counter_1 < tot_canv && abs(channel_info[sc_channel].peak_pos-channel_info[looK_channel].peak_pos) < 60)
            {
                TGraph *gr = new TGraph(waveform_size,samples,wf);
                graph_like_ivashkin_wants_it(gr,"Time [a.u.]","ADC channels", Form("center_waveform_%i_ch_%i",global_counter_1,channel_number),2);
                gr->GetXaxis()->SetRangeUser(200,350);
                gr->Draw("AL");
                global_counter_1++;
                flag++;
            }
*/
            if (global_counter_2 < tot_canv  && abs(channel_info[sc_channel].peak_pos-channel_info[looK_channel].peak_pos) < 400 
            && abs(channel_info[sc_channel].peak_pos-channel_info[looK_channel].peak_pos) > 192)
            {
                TGraph *gr = new TGraph(waveform_size,samples,wf);
                graph_like_ivashkin_wants_it(gr,"Time [a.u.]","ADC channels", Form("low_ampl_waveform_%i_ch_%i",global_counter_2,channel_number),2);
                gr->GetXaxis()->SetRangeUser(200,350);

                gr->Draw("AL");

                global_counter_2++;
                flag++;
            }
/*
            if (global_counter_3 < tot_canv  && abs(channel_info[sc_channel].peak_pos-channel_info[looK_channel].peak_pos) < 2000 
            && abs(channel_info[sc_channel].peak_pos-channel_info[looK_channel].peak_pos) > 960)
            {
                TGraph *gr = new TGraph(waveform_size,samples,wf);
                graph_like_ivashkin_wants_it(gr,"Time [a.u.]","ADC channels", Form("pedestal_waveworm_%i_ch_%i",global_counter_3,channel_number),2);
                gr->GetXaxis()->SetRangeUser(200,350);

                gr->Draw("AL");

                global_counter_3++;
                flag++;
            }
*/             
            if (flag > 0 && global_counter_2
             //+ global_counter_1 + global_counter_3 
             < 3*tot_canv) canv->SaveAs(source_path+"waveforms.pdf(","pdf");

            if (flag > 0 && global_counter_2 
            //+ global_counter_1 + global_counter_3 
            == 3*tot_canv ) canv->SaveAs(source_path+"waveforms.pdf)","pdf");
            delete canv;
            }
            }
            else
            {entryNum = total_entries;}
            }
        }

}