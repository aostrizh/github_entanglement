#ifndef CHANNEL_ENTRY_H
#define CHANNEL_ENTRY_H
#include<TTree.h>
#include <TNtuple.h>
#include <map>

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


struct mini_tree
{

    Float_t EdepScat0;
    Float_t EdepScat1;
    Float_t EdepDet0;
    Float_t EdepDet1;    
    Short_t DetNum0;
    Short_t DetNum1;
    Float_t EdepWeak;

        std::map<TString, Float_t*> branchfloat = 
        {
            {"EdepScat0", &EdepScat0},
            {"EdepScat1", &EdepScat1},
            {"EdepDet0", &EdepDet0},
            {"EdepDet1", &EdepDet1},     
            {"EdepWeak", &EdepWeak}        
        }; 
        std::map<TString, Short_t*> branchshort = 
        {
            {"DetNum0", &DetNum0},
            {"DetNum1", &DetNum1}
        }; 


    Int_t CreateBranches(TTree* mini_tree)
    {
        for (auto& val:branchfloat)
        {
            mini_tree->Branch(val.first,val.second,val.first+"\\F");
        }
        for (auto& val:branchshort)
        {
            mini_tree->Branch(val.first,val.second,val.first+"\\S");
        }

        return 1;
    }

    

    Int_t Initialize()
    {
        for (auto& val:branchfloat)        
        {
            
            val.second = 0;
        }
        for (auto& val:branchshort)
        {
            val.second = 0;

        }
        return 1;
    }

};

#endif
