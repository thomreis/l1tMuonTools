# l1tMuonTools
L1T muon tools for efficiency and rate studies from L1TNuples

##Usage examples:

###T&P efficiencies with batch system:
```
python create_batch_job.py -s muonTagAndProbe.py -p muonTagAndProbe -l input_file_list.txt -o ugmt_tandp_eff_histos --json json_file -w work_dir -j 10 -q 8nh --split_by_file
cd work_dir
source submit.sh
hadd ./out/ugmt_tandp_eff_histos.root ugmt_tandp_eff_histos_*.root
```
Produce efficiency plots
```
python plotTPEff.py -f work_dir/out/ugmt_tandp_eff_histos.root plotTPEff --eff --2d --emul
python plotTPEff.py -f work_dir/out/ugmt_tandp_eff_histos.root plotTPEff --eff --2d --run 123456
python plotTPEff.py -f work_dir/out/ugmt_tandp_eff_histos.root plotTPEff --fname2 work_dir2/out/ugmt_tandp_eff_histos.root --leg-txt1 text1 --leg-txt2 text2

```

###Rate histograms:
```
python makeSimpleRateHistos.py -f l1ntuple.root makeRateHistos
python plotSimpleRates.py -f ugmt_rate_histos.root plotRates
```
