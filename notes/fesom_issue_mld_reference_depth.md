# Draft FESOM issue: MLD2/MLD3 reference depth

Ready to file at https://github.com/FESOM/fesom2/issues/new
Title: MLD2/MLD3 reference the top model level, not 10 m, so they are not directly comparable with WOA18 mixed-layer depth

Evidence produced by `scripts/analysis/mld_definition_audit.py` (ARM=11X, YEAR=1389).

---

Checked against `main` @ a62f1806.

## What the code does

`src/oce_ale_pressure_bv.F90:462` computes MLD2 as the depth at which potential density exceeds its value at `nzmin`, the top model level:

```fortran
if ((rhopot(nz)-rhopot(nzmin) > sigma_theta_crit) .and. flag2) then
```

with `sigma_theta_crit = 0.125` kg/m3 (line 219). MLD3 is the same with `sigma_theta_crit_cmor = 0.03` (line 470).

## What the observational product does

The WOA18 mixed-layer climatology, the obvious reference to evaluate against, states its definition in the file's own `summary` attribute:

> Mixed layer is calculated for each profile by estimating the depth for which the potential density at 10m (reference depth) increases by 0.125 kg*m-3.

The threshold and the variable agree with MLD2 exactly. The reference depth does not: WOA uses 10 m, FESOM uses the top model level. On a mesh with a 5 m top layer that is 2.5 m. de Boyer Montegut et al. (2004), the other climatology in common use, also references 10 m, precisely to avoid the near-surface stratification a shallower reference picks up.

A related detail at line 426: MLD2 and MLD3 are initialised to `Z_3d_n(nzmin+1,node)`, so they can never return a value shallower than the second level (7.5 m here), whereas a 10 m-referenced MLD cannot return less than 10 m by construction.

## What it costs

Recomputed MLD from the model's own monthly T/S under both conventions. One coupled year, TCO95L91/CORE3, 220509 nodes, 47 levels, 5 m top layer, cavity nodes excluded, area-weighted:

| band | top-level reference | 10 m reference | difference |
|---|---|---|---|
| 90-60S | 365.5 m | 366.5 m | +1.05 m (0.3%) |
| 60-45S | 152.3 m | 152.8 m | +0.45 m (0.3%) |
| 45-30S | 104.0 m | 105.9 m | +1.95 m (1.9%) |
| 30S-30N | 41.9 m | 43.9 m | +1.98 m (4.7%) |
| 30-45N | 81.0 m | 83.4 m | +2.39 m (2.9%) |
| 45-60N | 137.4 m | 140.6 m | +3.15 m (2.3%) |
| 60-90N | 106.4 m | 109.9 m | +3.57 m (3.4%) |
| global | 91.1 m | 93.0 m | +1.93 m (2.1%) |

The offline recomputation under FESOM's own convention reproduces the model's MLD2 to a median of 1.2 m, the residual being EOS-80 against Jackett-McDougall, so the convention above is read correctly.

The effect is small in absolute terms and largest in relative terms where the mixed layer is shallow. It is also one-signed: MLD2 is biased shallow against any 10 m-referenced product everywhere, so it does not average out of a model-minus-observation map.

## Scope

This is a diagnostic comparability question, not a model physics one. MLD2 and MLD3 are written out and nothing else reads them. The field that feeds the GM parameterisation is `MLD1_ind` (`src/oce_fer_gm.F90:403`), which is unaffected.

## Two questions

1. Would you take a patch making the reference depth for MLD2/MLD3 a namelist parameter defaulting to 10 m, with the search starting at the first level below it? That makes both directly comparable with WOA18 and de Boyer Montegut without changing anything else.

2. MLD3 is labelled "Griffies" and the comment at line 422 says it exists "to be CMOR compliant". The CMIP6 `Omon` entry for `mlotst` says only "Sigma T is potential density referenced to ocean surface" and does not pin the reference depth or the threshold, so I cannot tell from the tables whether the current choice is intended to be protocol-conforming. If MLD3 is meant to be published as `mlotst`, the reference depth is worth stating in the field metadata either way.

## Unrelated detail noticed in the same block

Lines 463 and 471 use `rhopot(1)` in the interpolation while the test on lines 462 and 470 uses `rhopot(nzmin)`:

```fortran
MLD2(node)=MLD2(node)+(Z_3d_n(nz,node)-MLD2(node))/(rhopot(nz)-rhopot(nz-1)+1.e-20)*(rhopot(1)+sigma_theta_crit-rhopot(nz-1))
```

These are the same thing for open ocean, where `nzmin = 1`. Under ice-shelf cavities `nzmin > 1` and `rhopot` is only filled from `nzmin` downward (line 285), so `rhopot(1)` is still 0 from the initialisation at line 278. The last factor then evaluates to roughly `-rhopot(nz-1)`, about -1027 instead of a number of order 0.1.

In our output this does not surface, because cavity nodes come out masked. It looks like it should be `rhopot(nzmin)` regardless. Happy to send that as a separate one-line PR if you agree.
