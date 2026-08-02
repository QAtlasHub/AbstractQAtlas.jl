# Design: `@bound` — bounds declared as statements, with roles

Status: **shipped**. `@bound` replaces `@inequality`, which is removed. All 22
bound declarations are migrated. The generated behaviour is unchanged —
`residual` is still the `≥ 0` slack, `check` is still `slack ≥ −atol`, `solve`
still returns the saturation value — so this is an *expressiveness* change, not
a semantics change.

## 1. What `@inequality` could not say

```julia
@inequality :quantum LiebRobinsonBound(v, v_LR::LiebRobinsonVelocity) = v_LR - v
```

The statement being made is `v ≤ v_LR`. Nothing in the declaration says so. The
direction lives only in **the sign of an expression**, and the roles — which
slot is the subject, which one bounds it — are not recorded anywhere. Three
consequences:

1. **"What bounds X?" was not answerable.** `relations_constraining(q)` is
   role-blind: it returns `LiebRobinsonBound` for `LiebRobinsonVelocity`, the
   quantity that *does* the bounding, exactly as it would for the quantity being
   bounded. There was no query for the asymmetry.
2. **A consumer's `direction` field drifts.** Downstream registries carry
   `direction = :upper` beside a bound edge, maintained independently of the
   inequality's orientation. Two encodings of one fact, with no check between
   them.
3. **Writing the slack is a chance to get the sign backwards**, and a
   sign-flipped slack is a bound that passes on every input — the failure mode
   is silence, not an error.

## 2. The three forms

```julia
# comparison — the parameter list IS the statement
@bound :thermodynamic SpecificHeatPositivity(Cv::SpecificHeat >= 0)
@bound :quantum LiebRobinsonBound(v <= v_LR::LiebRobinsonVelocity)

# statement — slots declared normally, body is a comparison over them
@bound :entanglement Subadditivity(S_A, S_B, S_AB) = S_AB <= S_A + S_B
@bound :quantum MandelstamTammBound(τ, ΔE) = τ >= π / (2 * ΔE)

# bare slack — absorbs @inequality verbatim; no direction claimed
@bound :test SomeBound(x, y) = x - y
```

**The subject is written first.** `C ≥ 0`, `v ≤ v_LR`, `S_AB ≤ S_A + S_B` all
read as sentences about their subject, and a constant on the left is a
declaration error telling you to flip the operator. This is what makes
`bounded_slot` unambiguous rather than a guess.

Strict `<` / `>` are **rejected**. The criterion is `slack ≥ −atol`, so a strict
declaration would be checked non-strictly — a silent weakening. State the
non-strict form.

## 3. Recorded roles

| trait | comparison form | statement form | bare slack | equality |
|---|---|---|---|---|
| `bounded_slot` | the left slot | the left slot *if it is one variable* | `nothing` | `nothing` |
| `bounding_slot` | the right slot | the right slot *if it is one variable* | `nothing` | `nothing` |
| `bounding_constant` | the right literal | the right literal | `nothing` | `nothing` |
| `bound_direction` | `:upper` / `:lower` | `:upper` / `:lower` | `:slack` | `nothing` |

The statement form records **only what it can**: `StrongSubadditivity`'s
`S_ABC + S_B ≤ S_AB + S_BC` has an expression on both sides, so both slots are
`nothing`. That is honest — there is no single subject — and better than
nominating one. `:slack` and `nothing` are distinct: the first means "a bound
that stated no direction", the second means "not a bound".

`bounds_on(q)` is the reverse index, matching on the *bounded* slot's identity
type. It is deliberately narrower than `relations_constraining`:
`bounds_on(LiebRobinsonVelocity)` is empty because that quantity does the
bounding, while `relations_constraining(LiebRobinsonVelocity)` still finds the
relation. Family slots resolve through it (§8a), so
`bounds_on(Susceptibility{(:z,:z)})` finds the family-declared
`SusceptibilityPositivity`.

## 4. The extension win: a bound may outrun its vocabulary

In the comparison form the **bounded side may be an untyped slot** while the
bounding side carries a quantity type:

```julia
@bound :quantum LiebRobinsonBound(v <= v_LR::LiebRobinsonVelocity)
```

`v` is any independently measured information velocity, and no quantity names it
— the atlases hold the bound, not the measurement. The same shape covers the
Tsirelson/Mermin/MSS/Bekenstein family, where the bounding VALUE is fetchable
and the bounded observable has no type in either package: the bound can be
stated, and its value fetched, before a quantity exists for what it bounds.

## 5. What deliberately stays out

**Saturation is not expressible, by decision.** That a particular model
saturates a bound (TFIM on the Lieb–Robinson velocity) is model-specific and
belongs on the consumer's registry row as `scheme = :saturating_bound`, not on
a universal relation. `solve` already returns the saturation *value*; asserting
that a model *attains* it is an atlas claim.

## 6. The silent-failure guard

Statement-form bodies reach a macro wrapped in a `:block` carrying a
`LineNumberNode`. Comparison detection must unwrap it — an un-unwrapped
`τ >= π/(2ΔE)` reads as an ordinary slack expression, the kernel then *evaluates*
the comparison instead of measuring it and returns `true`, and `true ≥ −atol`
holds for every input. A bound that can never fail, with no error anywhere. This
was a real bug during the migration, caught by probe rather than by construction.

Two things now make it non-silent: the unwrap, and a runtime guard that rejects a
`Bool` slack at `residual` with an explanation. `test_bound_macro.jl` sweeps the
whole registry asserting every slack is a `Number` — the structural version of the
per-relation saturation cases, which cover today's bounds but not tomorrow's.
