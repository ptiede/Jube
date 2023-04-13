# Jube API



## Contents

```@contents
Pages = ["api.md"]
```

## Index

```@index
Pages = ["api.md"]
```

## Metrics

### Abstract Metric

```@docs
Jube.metric_uu
Jube.metric_dd
```

### Kerr

```@docs
Jube.metric_uu
Jube.metric_dd
Jube.horizon
Jube.rs
Jube.η    
Jube.λ
```


### Combinators

```@docs
Base.:+(::Comrade.AbstractModel, ::Comrade.AbstractModel)
Comrade.added
Comrade.convolved
Comrade.components
Comrade.smoothed
Comrade.CompositeModel
Comrade.AddModel
Comrade.ConvolvedModel
```

### Geometric and Image Models

```@docs
Comrade.GeometricModel
Comrade.ConcordanceCrescent
Comrade.Crescent
Comrade.Disk
Comrade.SlashedDisk
Comrade.ExtendedRing
Comrade.Gaussian
Comrade.MRing
Comrade.Ring
Comrade.ParabolicSegment
Comrade.ContinuousImage
Comrade.ZeroModel
```

### Image Pulses
```@docs
Comrade.Pulse
Comrade.DeltaPulse
Comrade.BSplinePulse
Comrade.RaisedCosinePulse
Comrade.BicubicPulse
Comrade.Butterworth
```


### Model Image (non analytic FFT)

```@docs
Comrade.create_cache
Comrade.update_cache
Comrade.modelimage
Comrade.uviterator
Comrade.fouriermap
Comrade.ModelImage
Comrade.DFTAlg
Comrade.DFTAlg(::Comrade.EHTObservation)
Comrade.DFTAlg(::AbstractArray, ::AbstractArray)
Comrade.DFTAlg(::Comrade.ArrayConfiguration)
Comrade.FFTAlg
Comrade.FFTCache
Comrade.NFFTAlg
Comrade.NUFTCache
Comrade.ObservedNUFT
```


### Modifiers

```@docs
Comrade.basemodel
Comrade.unmodified
Comrade.renormed
Comrade.rotated
Comrade.posangle
Comrade.shifted
Comrade.stretched
Comrade.ModifiedModel
Comrade.AbstractModifier
Comrade.Stretch
Comrade.Shift
Comrade.Rotate
Comrade.Renormalize
```

### Polarized Models

```@docs
ComradeBase.mbreve
ComradeBase.evpa(pimg::ComradeBase.AbstractPolarizedModel, p)
Comrade.PolarizedModel
```


### Model Evaluation

For more docstrings on how to evaluate models see [ComradeBase API](@ref).

```@docs
Comrade.amplitude
Comrade.amplitudes
Comrade.bispectra
Comrade.bispectrum
Comrade.closure_phase
Comrade.closure_phases
Comrade.logclosure_amplitude
Comrade.logclosure_amplitudes
Comrade.visibility
```

## Data Types


```@docs
Comrade.amplitude(::Comrade.EHTVisibilityDatum)
Comrade.amplitude(::Comrade.EHTVisibilityAmplitudeDatum)
Comrade.baselines
Comrade.arrayconfig
Comrade.closure_phase(::Comrade.EHTVisibilityDatum, ::Comrade.EHTVisibilityDatum, ::Comrade.EHTVisibilityDatum)
Comrade.getdata
Comrade.getuv
Comrade.getuvtimefreq
Comrade.scantable
Comrade.stations
Comrade.uvpositions
Comrade.ArrayConfiguration
Comrade.ClosureConfig
Comrade.AbstractInterferometryDatum
Comrade.ArrayBaselineDatum
Comrade.EHTObservation
Comrade.EHTArrayConfiguration
Comrade.EHTCoherencyDatum
Comrade.EHTClosurePhaseDatum
Comrade.EHTLogClosureAmplitudeDatum
Comrade.EHTVisibilityDatum
Comrade.EHTVisibilityAmplitudeDatum
Comrade.Scan
Comrade.ScanTable
```

## eht-imaging interface

```@docs
Comrade.extract_coherency
Comrade.extract_amp
Comrade.extract_cphase
Comrade.extract_lcamp
Comrade.extract_vis
Comrade.load_ehtim_uvfits
Comrade.load_ehtim
Comrade.scan_average
```

## Bayesian Tools

### Posterior Constructions

```@docs
Comrade.ascube
Comrade.asflat
Comrade.flatten
Comrade.inverse
Comrade.prior_sample
Comrade.likelihood
Comrade.sample(::Posterior)
Comrade.transform
Comrade.MultiRadioLikelihood
Comrade.Posterior
Comrade.TransformedPosterior
Comrade.RadioLikelihood
Comrade.IsFlat
Comrade.IsCube
```

## Misc

```@docs
Comrade.μas2rad
Comrade.rad2μas
Comrade.load
Comrade.save
Comrade.NonAnalyticTest
Comrade.station_tuple
```

## Internal (Not Public API)
```@docs
Comrade.scale_uv
Comrade.scale_image
Comrade.transform_uv
Comrade.transform_image
Comrade.ThreadedModel
Comrade.extract_FRs
```
