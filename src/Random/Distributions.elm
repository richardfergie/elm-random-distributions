module Random.Distributions (normal) where

{-| This library provides standard non-uniform random sampling methods for the
core Random library.

# Distributions implemented
@docs normal

-}

import Random

{-| The error function
Approximation with a maximal error of 1.2*10^-7.

Directly from wikipedia:
https://en.wikipedia.org/wiki/Error_function#Numerical_approximation

\begin{align}
\tau = {} & t\cdot\exp\left(-x^2-1.26551223+1.00002368 t+0.37409196 t^2+0.09678418 t^3\right.\\
& \left.{}-0.18628806 t^4+0.27886807 t^5-1.13520398 t^6+1.48851587\cdot t^7\right. \\
& \left.{}-0.82215223 t^8+0.17087277 t^9\right)
\end{align}
-}
erf x =
  let
    t = 1 / (1 + 0.5 * abs x)
    exponent =
      -x^2
      - 1.26551223
      + 1.00002368 * t
      + 0.37409196 * t^2
      + 0.09678418 * t^3
      - 0.18628806 * t^4
      + 0.27886807 * t^5
      - 1.13520398 * t^6
      + 1.48851587 * t^7
      - 0.82215223 * t^8
      + 0.17087277 * t^9
    tau = t * e ^ exponent
    y =
      if x >= 0
        then 1 - tau
        else tau + 1
  in
    clamp 0 1 y

{-| The complimentary error function
Approximation with a maximal error of 1.2*10^-7.
-}
erfc x = clamp 0 1 (1 - erf x)

{-| The natural logarithm
-}
ln x = logBase e x

probability : Random.Generator Float
probability =
  Random.float 0 1

{-| A `random' generator that always returns the given value.
-}
identityGenerator : Float -> Random.Generator Float
identityGenerator x =
  Random.map (always x) Random.bool

-- probabilities : Generator (Float, Float)
-- probabilities = pair probability probability

normalDensity mu sigma x =
  let
    factor = 1 / (sigma * sqrt (2*pi))
    exponent = -(x-mu)^2 / (2*sigma^2)
  in
    factor * e ^ exponent

normalDensityInverse mu sigma y =
  let
    factor = 1 / (sigma * sqrt (2*pi))
  in
    mu + sqrt (-2 * sigma^2 * ln (y / factor))


tableSize = 256
normalZigguratTable = [0..tableSize]

{-| Implement the Ziggurat algorithm.

https://en.wikipedia.org/wiki/Ziggurat_algorithm

-}
ziggurat density = Random.float 0 1
--   let n =

{-| Fallback algorithm for the tail of a normal distribution

From wikipedia: https://en.wikipedia.org/wiki/Ziggurat_algorithm

For a normal distribution, Marsaglia suggests a compact algorithm:
  1.  Let x = −ln(U1)/x1.
  2.  Let y = −ln(U2).
  3.  If 2y > x^2, return x + x1.
  4.  Otherwise, go back to step 1.
-}
zigguratNormalTail x1 =
  let
    u1u2gen = Random.pair probability probability
    fallback (u1, u2) =
      let
        x = -(ln u1)/x1
        y = -(ln u2)
      in
        if 2*y > x^2
          then identityGenerator (x + x1)
          else zigguratNormalTail x1
  in
    u1u2gen `Random.andThen` fallback

{-| Generate a standard normal distribution using the Ziggurat algorithm.

https://en.wikipedia.org/wiki/Ziggurat_algorithm

-}
normal : Random.Generator Float
-- normal = ziggurat <| normalDensity 0 1
normal = zigguratNormalTail 3.5
