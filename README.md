# l1tMuonTools
L1T muon tools for efficiency, rate and other studies from L1TNuples

## Efficiencies
Numerator and denominator histograms for tag and probe efficiencies are calculated with the `muonTagAndProbe.py` script. L1 ntuples from a SingleMuon dataset (eventually a ZMu skim) are taken as the input.
```
python muonTagAndProbe.py -l input_l1ntuple_file_list.txt muonTagAndProbe --json good_ls_json.txt --outname ugmt_tandp_eff_histos.root --era 2017pp --use-l1-extra-coord
```
With the `-l` option a text file with the paths to the input files can be used and with the `-f` option a single L1 ntuple input file can be specified.
To run on emulated muons add the `--emul` option. With the `--run` option a list of runs to be analysed can be selected.

### Using the batch system:
To run over many input files the task can be divided and sent to the lxbatch system.
```
python create_batch_job.py -s muonTagAndProbe.py -p muonTagAndProbe -l input_l1ntuple_file_list.txt -w work_dir --cmd-line-args " --json good_ls_json.txt --outname ugmt_tandp_eff_histos.root --era 2017pp --emul" -njobs 10 -queue 8nh --split_by_file --submit
```
Once the jobs are finished combine them with `hadd`.
```
hadd ./work_dir/out/ugmt_tandp_eff_histos.root ugmt_tandp_eff_histos_*.root
```

### Calculating the integrated luminosity
The integrated luminosity can be calculated with the `brilcalc` tool and the processedLumis json from `crab report`.
If the intersection of two LS json files is needed (For example of the processedLumis json with the golden json) the `compareJSON.py` script available in CMSSW can be used with the `--and` option.

### Produce efficiency plots
Once the ROOT file with the numerator and denominator histograms is generated the efficiency plots can be produced with the `plotTPEff.py` script.
```
python plotTPEff.py -f ugmt_tandp_eff_histos.root plotTPEff --year 2016 --lumi="33.0 fb^{-1}" --eff --2d --qualcomp --delta --control
```
* `--eff` produces the efficiency plots. To run on emulated muons used the `--emul` option.
* `--2d` produces 2D plots with RECO vs. L1 for several muon variables.
* `--qualcomp` produces efficiency plots for different minimal L1 muon qualities in one plot.
* `--delta` produces plots with the difference between RECO probe muon and matched L1 muon.
* `--control` produces control histograms for the probe and L1 muons.
* `--public` produces publication style plots and output files in png and pdf format.
* `--data-emul` produces efficiency comparison plots between data and emulator. To use this option the `muonTagAndProbe.py` script has to be run with and without the `--emul` option, and the two ROOT files have to be merged with `hadd`.
* `--upgrade-legacy` produces efficiency comparison plots between upgrade and legacy. To use this option the `muonTagAndProbe.py` script has to be run with and without the `--legacy` option, and the two ROOT files have to be merged with `hadd`.
* Efficiencies for one specific run can be plotted with the `--run` option.
* Comparing efficiencies from different ugmt_tandp_eff_histos.root files can be done with the `--fname2` option that specifies the second input file. This can also be used to compare two runs with the `--run2` option.

```
python plotTPEff.py -f ugmt_tandp_eff_histos.root plotTPEff --fname2 ugmt_tandp_eff_histos2.root --leg-txt1 text1 --leg-txt2 text2
python plotTPEff.py -f ugmt_tandp_eff_histos.root plotTPEff --fname2 ugmt_tandp_eff_histos.root --run 123456 --run2 123457 --leg-txt1 run123456 --leg-txt2 run123457
```

## Rates:
Rates are calculated from ZeroBias samples with the `makeSimpleRateHistos.py` script. The same options `-l` and `-f` as for the tag and probe tool are available.
```
python makeSimpleRateHistos.py -f l1ntuple.root makeRateHistos --json good_ls_json.txt --use-l1-extra-coord
```
* `--emul` uses emulated muons
* `--use-l1-extra-coord` uses the eta and phi coordinates extrapolated to the vertex by the uGMT

### Using the batch system:
The batch submission script `create_batch_job.py` can be used as well to split the task in more jobs and send them to lxbatch.
```
python create_batch_job.py -s makeSimpleRateHistos.py -p makeRateHistos -l input_l1ntuple_file_list.txt -w work_dir --cmd-line-args " --json good_ls_json.txt --use-l1-extra-coord" -njobs 10 -queue 8nh --split_by_file --submit
```

### Produce rate plots:
To produce the rate plots with the correct rate in kHz the number of colliding bunches in CMS is needed as an input. This is the third number in the LHC filling scheme (e.g. fill 5966: 25ns_2556b_2544_2215_2332_144bpi_20inj has 2544 bunches colliding in CMS).
```
python plotSimpleRates.py -f ugmt_rate_histos.root plotRates --bunches 2544
```
* `--scale` defines an additional scale factor
* `--scale-to-bunches` linearly extrapolates the rates from the bunches given with `--bunches` to the new number of bunches.
* To compare histograms from different files or runs the `--fname2` option can be used.

In order to get the rate for a MC sample the number of bunches corresponding to a desired instantaneous luminosity and PU needs to be calculated by the tool. Therefore, the options `--instlumi`, `--pu`, and `--xsect` need to be used, where the last one defines the total cross section in mb.

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
  python create_batch_job.py -l files_SingleMu_Pt1To1000_FlatRandomOneOverPt_bmtfetafix.txt -q 8nh -j 25 --split_by_file -s muonExtrapolation.py -p muonExtrapolation -w muon_extrapolation_histos_SingleMu_pt1To1000_i92p9_5etabits_posside_bmtf --cmd-line-args " -o muon_extrapolation_histos --emul --eta-bits 5 --tf-type b --pos-side" --submit

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
