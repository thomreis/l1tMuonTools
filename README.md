# l1tMuonTools
L1T muon tools for efficiency, rate and other studies from L1TNuples

## T&P efficiencies with batch system:
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

## Rate histograms:
```
python makeSimpleRateHistos.py -f l1ntuple.root makeRateHistos
python plotSimpleRates.py -f ugmt_rate_histos.root plotRates
```

## Muon coordinates at the vertex
The uGMT can project the muon eta and phi coordinates, measured usually at the 2nd muon station, to the vertex by means of LUTs.

Workflow to produce these LUTs

#### Take a L1Ntuple from a SingleMuon gun MC sample with L1 muon and GEN information
#### Produce the histograms with the delta beween L1 measurement and GEN
  * The histograms can be produced for each track finder separately with the `--tftype {b|o|e}` option (BMTF (b), OMTF (o), EMTF (e)). Combinations are possible as well.
  * The emulated branches of the ntuples are used with the `--emul` option.
  * The `--eta-bits` option specifies the number of eta bits that are going to be used in the LUT and, thus, for how many eta ranges histograms have to be produced.
  * Since the eta projection depends on the detector side one can use the `--pos-eta` or `--neg-eta` options to select only one detector side.
  ```
  python muonExtrapolation.py -f /path/to/L1Ntuple.root muonExtrapolation --emul --eta-bits 5
  ```
  - Running on the lxbatch system to speed things up
  ```
  python create_batch_job.py -l files_SingleMu_Pt1To1000_FlatRandomOneOverPt_bmtfetafix.txt -q 8nh -j 25 --split_by_file -s muonExtrapolation.py -p muonExtrapolation -w muon_extrapolation_histos_SingleMu_pt1To1000_i92p9_5etabits_posside_bmtf -o muon_extrapolation_histos --emul --eta-bits 5 --tf-type b --pos-side
  cd muon_extrapolation_histos_SingleMu_pt1To1000_i92p9_5etabits_posside_bmtf/
  source submit.sh
  ```
  - The resulting histograms can be plotted
  ```
  python plotExtrapolation.py -f muon_extrapolation_histos.root plotExtrapolation --delta --2d -i
  ```

#### Produce the LUTs for the different track finders separately
  ```
  python generateExtrapolationLut.py -f muon_extrapolation_histos.root generateExtrapolationLut --coordinate {eta|phi} --pt-bits 7 --eta-bits 5 --out-bits 4 --out-shift 2
  ```
  * Input widths: `--pt-bits` and `--eta-bits`
  * Output width: `--out-bits`
  * Bit shift of the output (essentially a factor 2^shift): `--out-shift`

  * The LUTs are kept in `/L1Trigger/L1TMuon/data/microgmt_luts/`

#### Use the new LUTs to generate L1Ntuples from a SingleMuon MC sample
#### Make and plot new histograms with the coordinates at the vertex
  * Using the new coordinates at the vertex with `--extrapolated`.
  ```
  python plotExtrapolation.py -f muon_extrapolation_histos.root plotExtrapolation --delta --2d --extrapolated -i
  ```
