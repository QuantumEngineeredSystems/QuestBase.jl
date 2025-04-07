```@meta
CollapsedDocStrings = true
```

# API

```@contents
Pages = ["API.md"]
Depth = 2:3
```

## Classical harmonic systems

### Classical equations of motion in lab frame

```@docs
QuestBase.DifferentialEquation
QuestBase.rearrange_standard!(::QuestBase.DifferentialEquation)
QuestBase.rearrange!
QuestBase.is_rearranged_standard
QuestBase.get_equations
QuestBase.get_independent_variables(::QuestBase.DifferentialEquation)
QuestBase.get_variables(::QuestBase.DifferentialEquation)
```

### Harmonics

```@docs
QuestBase.HarmonicVariable
QuestBase.add_harmonic!
```

### Effective stroboscopic equations of motion in a rotating frame

```@docs
QuestBase.HarmonicEquation
QuestBase.get_independent_variables(::QuestBase.HarmonicEquation)
QuestBase.get_variables(::QuestBase.HarmonicEquation)
QuestBase.rearrange_standard(::QuestBase.HarmonicEquation)
QuestBase.rearrange_standard!(::QuestBase.HarmonicEquation)
QuestBase.rearrange
QuestBase.is_rearranged
```

## Steady state methods

```@docs
QuestBase.HarmonicBalanceMethod
```

## Symbolic Utilities

QuestBase contains a number of symbolic utilities to help with the symbolic manipulation of the equations of motion. These are function on the top of the Symbolics.jl package and are considered **non-public**.

```@docs
QuestBase.d
```

### Exponentials

```@docs
QuestBase.simplify_exp_products
QuestBase.expand_exp_power
```

### Trigonometrics

```@docs
QuestBase.trig_reduce
QuestBase.is_trig
QuestBase.exp_to_trig
QuestBase.trig_to_exp
```