#ifndef CHANNEL_ENTRY_H
#define CHANNEL_ENTRY_H
#include<TTree.h>

const int MAX_N_SAMPLES = 2048;

/*
struct ChannelEntry {
    qint32 integral_in_gate;
    qint16 zero_level;
    qint16 peak_pos;
    qint16 amp;
    qint16 wf_size;
    qint16 wf[MAX_N_SAMPLES];

};
*/
struct ChannelEntry {
    Int_t integral_in_gate;
    Short_t zero_level;
    Short_t peak_pos;
    Short_t amp;
    Short_t wf_size;
    Short_t wf[MAX_N_SAMPLES];

    static TString GetChName(Int_t channel_num)
    {
	return TString::Format("channel_%i", channel_num);
    }

    Int_t SetBranch(TTree *tree, Int_t channel_num)
    {
	return tree->SetBranchAddress(GetChName(channel_num).Data(), this);
    }
};

struct short_ChannelEntry
{
    Int_t integral_in_gate;
    Short_t peak_pos;
    UShort_t amp; 

    static TString GetChName(Int_t channel_num)
    {
	return TString::Format("channel_%i", channel_num);
    }

    TBranch* CreateBranch(TTree *tree, Int_t channel_num)
    {
	    return tree->Branch(GetChName(channel_num).Data(), this, "integral_in_gate/I:peak_pos/S:amp/s");
    }
        Int_t SetBranch(TTree *tree, Int_t channel_num)
    {
	return tree->SetBranchAddress(GetChName(channel_num).Data(), this);
    }
};

    struct short_energy_ChannelEntry
{
    Float_t integral_in_gate;
    Short_t peak_pos;
    UShort_t amp; 

    static TString GetChName(Int_t channel_num)
    {
	return TString::Format("channel_%i", channel_num);
    }

    TBranch* CreateBranch(TTree *tree, Int_t channel_num)
    {
	    return tree->Branch(GetChName(channel_num).Data(), this, "integral_in_gate/F:peak_pos/S:amp/s");
    }
        Int_t SetBranch(TTree *tree, Int_t channel_num)
    {
	return tree->SetBranchAddress(GetChName(channel_num).Data(), this);
    }
    
};

#endif
