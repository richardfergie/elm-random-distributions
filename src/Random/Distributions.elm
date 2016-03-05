module Random.Distributions
  ( erf
  , erfc
  , erfinv
  , normal
  , normalCumulative
  , normalDensity
  , normalDensityInverse
  , normalQuantile
  , probit
  , quantile
  -- , simpleNormal
  , ziggurat
  , zigguratTables
  , zigguratX1
  ) where

{-| This library provides non-uniform distributions for the core Random library.

# Generators
@docs normal

# Distribution functions
@docs normalCumulative
@docs normalDensity
@docs normalDensityInverse
@docs normalQuantile
@docs erf
@docs erfc
@docs erfinv
@docs probit

# Other functions
@docs quantile

## Ziggurat algorithm

Helper functions implementing the [Ziggurat
algorithm](https://en.wikipedia.org/wiki/Ziggurat_algorithm).

@docs ziggurat
@docs zigguratTables
@docs zigguratX1

-}

import Array
import Random
import String

{-| Return elements of the given list at the given indexes.
Assumes indexes are sorted.
Skips out of range indexes.
-}
getIndexes : List a -> List number -> List a
getIndexes list indexes =
  let
    do index xs is =
      case (List.head xs, List.head is) of
        (Just x, Just i) ->
          let
            newXs = List.drop 1 xs
            newIndex = index + 1
          in
            if i == index
              then
                x :: do newIndex newXs (List.drop 1 is)
              else
                do newIndex newXs is
        a ->
            []
  in
    do 0 list indexes

{-| Produces sample quantiles of the xs corresponding to the given probabilities.

    quantile samples probs

Based on the [corresponding R function](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/quantile.html).

`probs` is assumed to be a sorted list of probabilities.
-}
quantile : List Float -> List Float -> List Float
quantile samples probs =
  let
    numSamples = List.length samples
    fracProbIndex p = p * (toFloat numSamples)
    fracIndices = List.map fracProbIndex probs

    maxIndex = numSamples - 1
    lowerIndices = List.map (clamp 0 maxIndex << floor) fracIndices
    upperIndices = List.map (clamp 0 maxIndex << ceiling) fracIndices

    sorted = List.sort samples
    lowerQuantiles = getIndexes sorted lowerIndices
    upperQuantiles = getIndexes sorted upperIndices

    fractions = List.map (\f -> f - toFloat (truncate f)) fracIndices

    interpolate f l u = (1-f)*l + f*u
  in
    List.map3 interpolate fractions lowerQuantiles upperQuantiles

-- -- Simple implementation of quantile without interpolation between samples.
-- quantile : List Float -> List Float -> List Float
-- quantile samples probs =
--   let
--     numSamples = List.length samples
--     probIndex p = floor <| p * (toFloat numSamples)
--     indexes = List.map probIndex probs
--     sorted = List.sort samples
--   in
--     getIndexes sorted indexes

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
erf : Float -> Float
erf x =
  let
    t = 1.0 / (1 + 0.5 * (abs x))
    exponent = (
      -(x^2)
      - 1.26551223
      + 1.00002368 * t
      + 0.37409196 * t^2
      + 0.09678418 * t^3
      - 0.18628806 * t^4
      + 0.27886807 * t^5
      - 1.13520398 * t^6
      + 1.48851587 * t^7
      - 0.82215223 * t^8
      + 0.17087277 * t^9)
    tau = t * e ^ exponent
    y =
      if x >= 0
        then 1 - tau
        else tau - 1
  in
    clamp -1 1 y

{-| The complimentary error function.
Approximation with a maximal error of 1.2*10^-7.
-}
erfc : Float -> Float
erfc x = clamp 0 2 (1 - erf x)

{-| The inverse of the error function.

Implementation [from wikipedia](https://en.wikipedia.org/wiki/Error_function#Approximation_with_elementary_functions)
-}
erfinv : Float -> Float
erfinv x =
  let
    sgn = if x > 0 then 1 else -1
    -- a = 0.147 -- error < 0.00012
    a = (8*(pi-3)) / (3*pi*(4 - pi)) -- error < 0.00035 very accurate near 0 and inf
    ln1_x2 = ln (1 - x^2)
    parens = 2/(pi*a) + ln1_x2 / 2
    sqrt1 = parens^2 - ln1_x2 / a
    terms = sqrt sqrt1 - parens
    result = sgn * sqrt terms
  in
    if x == 0
      then 0
      else result

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

{-| The probit function.

The quantile function for the standard normal distribution (i.e. the inverse of
the cumulative distribution function of the standard normal distribution).

Implemented using the inverse error function [as described on Wikipedia](https://en.wikipedia.org/wiki/Probit#Computation)
-}
probit : Float -> Float
probit p =
  (sqrt 2) * erfinv (2*p - 1)

{-| The quantile function for a normal distribution with the given mean and standard deviation

i.e. the inverse of the cumulative distribution function of the normal
distribution.

    q = normalQuantile mu sigma p

Implementation using the probit function [from Wikipedia](https://en.wikipedia.org/wiki/Normal_distribution#Quantile_function).
-}
normalQuantile : Float -> Float -> Float -> Float
normalQuantile mu sigma p =
  mu + sigma * probit p


{-| The cumulative distribution function of the standard normal distribution.

    y = standardNormalCdf x

-}
standardNormalCumulative : Float -> Float
standardNormalCumulative x =
  0.5 * (1 + erf (x / sqrt 2))

{-| The cumulative distribution function of a normal distribution.

Implemented using [the error function](https://en.wikipedia.org/wiki/Normal_distribution#Cumulative_distribution_function).

    y = normalCumulative mu sigma x

-}
normalCumulative : Float -> Float -> Float -> Float
normalCumulative mu sigma x =
  standardNormalCumulative <| (x - mu) / sigma

{-| The probability density function for a normal distribution.

    y = normalDensity mu sigma x

-}
normalDensity : Float -> Float -> Float -> Float
normalDensity mu sigma x =
  let
    factor = 1 / (sigma * sqrt (2*pi))
    exponent = -((x-mu)^2 / (2*sigma^2))
  in
    factor * e ^ exponent

{-| The inverse of the density function for a normal distribution.

    x = normalDensityInverse mu sigma y

-}
normalDensityInverse : Float -> Float -> Float -> Float
normalDensityInverse mu sigma y =
  let
    z = max 0 y
    factor = 1 / (sigma * sqrt (2*pi))
    scaledProb = max 0 <| z / factor
    sqrDeviance = max 0 <| (-2 * sigma^2 * ln (scaledProb))
  in
    mu + sqrt sqrDeviance

{-| Find x1 and A for a given table size, density function, and inverse density
function.

https://en.wikipedia.org/wiki/Ziggurat_algorithm#Finding_x1_and_A

    x1 = zigguratX1 n pdfFunc invPdfFunc cdfFunc

-}
zigguratX1 : Int -> (Float -> Float) -> (Float -> Float) -> (Float -> Float) -> Float
zigguratX1 n pdfFunc invPdfFunc cdfFunc =
  let
    f0 = pdfFunc 0
    areaDiffFunc x1 =
      let
        y1 = pdfFunc x1
        tailArea = 1 - cdfFunc x1
        baseLayerArea = x1*y1 + tailArea
        tables = zigguratTables n y1 baseLayerArea pdfFunc invPdfFunc
        (xn_1, yn_1) =
          case List.head <| List.drop n <| tables of -- get the last element
            Just pair -> pair
            Nothing -> Debug.crash "The list normal ziggurat tables was not of length n"
        topLayerArea = xn_1*(f0 - yn_1)
      in
        topLayerArea - baseLayerArea

    searchEps = 1e-12
    upperBound = 5

    -- Find the lowest x1 for which areaDiffFunc returns a valid result.
    -- Without this, the bisection search for x1 will fail.
    diffValid x1 =
      let
        diff = areaDiffFunc x1
      in
        if isInfinite diff || isNaN diff
          then -1 else 1
    lowerBound = case bisectionSearch diffValid searchEps 100 0 upperBound of
      Just v -> v + searchEps
      Nothing -> Debug.crash "Could not find a stable lower bound for x1 while generating the normal ziggurat tables."

    -- Perform a bisection search for x1.
    x1 =
      case bisectionSearch areaDiffFunc searchEps 100 lowerBound upperBound of
        Just v -> v
        Nothing -> Debug.crash "The bisectionSearch failed while generating the normal ziggurat tables."
  in
    x1

{-| Bisection method for root finding

https://en.wikipedia.org/wiki/Bisection_method
-}
bisectionSearch : (Float -> Float) -> Float -> Int -> Float -> Float -> Maybe Float
bisectionSearch f eps n a b =
  let
    sign x =
      if x > 0
        then 1
        else -1
    search n a b =
      let
        va = f a
        vb = f b
      in
        if n <= 0 || (sign va == sign vb)
          then Nothing
          else
            let
              c = (a + b) / 2
              vc = f c
            in
              if vc == 0 || (b - a) / 2 < eps
                then Just c
                else
                  if sign vc == sign va
                    then search (n-1) c b
                    else search (n-1) a c
  in
    if a < b
      then search n a b
      else search n b a


{-| Generate the ziggurat tables.

https://en.wikipedia.org/wiki/Ziggurat_algorithm#Generating_the_tables

    tables = zigguratTables numLayers y1 layerArea pFunc invPFunc
    (xs, ys) = List.unzip tables

-}
zigguratTables : Int -> Float -> Float -> (Float -> Float) -> (Float -> Float) -> List (Float, Float)
zigguratTables n y1 layerArea pFunc invPFunc =
  let
    x1 = invPFunc y1
    nextLayer (xi, yi) =
      let
        yi1 = yi + layerArea / xi
        xi1 = invPFunc yi1
      in (xi1, yi1)
    layerList = List.scanl (\_ x1y1 -> nextLayer x1y1) (x1, y1) [1..n]
    ficticiousX0Y0 = (layerArea/y1, 0)
  in
    ficticiousX0Y0 :: layerList



{-| Implement the [Ziggurat algorithm](https://en.wikipedia.org/wiki/Ziggurat_algorithm) for one-sided distributions.

    oneSidedNormalGenerator =
      let
        n = numLayers
        pdfFunc = normalDensity 0 1
        invPdfFunc = normalDensityInverse 0 1
        cdfFunc = standardNormalCumulative
        x1 = zigguratX1 n pdfFunc invPdfFunc cdfFunc
        y1 = pdfFunc x1
        listTables = zigguratTables n y1 layerArea pdfFunc invPdfFunc
        tables = Array.fromList listTables
        tailGen = zigguratNormalTail x1
      in
        ziggurat tables pdfFunc tailGen

-}
ziggurat : Array.Array (Float, Float) -> (Float -> Float) -> Random.Generator Float -> Random.Generator Float
ziggurat tables pFunc tailGen =
  let
    numLayers = Array.length tables
    layerGen = Random.int 0 (numLayers-2)
    layerU0U1gen = Random.map3 (,,) layerGen probability probability

    chooseLayer (i, u0, u1) =
      let
        (xi, yi) = case Array.get i tables of
          Just pair -> pair
          Nothing -> Debug.crash (String.append "Index i='" <| String.append (toString i) "'out of range")
        (xi1, yi1) = case Array.get (i+1) tables of
          Just pair -> pair
          Nothing -> Debug.crash (String.append "Index i+1='" <| String.append (toString (i+1)) "'out of range")
        x = u0*xi
      in
        if x < xi1
          then
            identityGenerator x
          else
            if i == 0
              then
                tailGen
              else
                let
                  y = yi + u1*(yi1 - yi)
                  fx = pFunc x
                in
                  if y < fx
                    then
                      identityGenerator x
                    else
                      layerU0U1gen `Random.andThen` chooseLayer
  in
    layerU0U1gen `Random.andThen` chooseLayer

{-| Fallback algorithm for the tail of a normal distribution.
-}
zigguratNormalTail : Float -> Random.Generator Float
zigguratNormalTail x1 =
  let
    p1 = standardNormalCumulative -x1
    fallback p = probit (p*p1)
  in
    Random.map fallback probability

{-| Fallback algorithm for the tail of a normal distribution

From wikipedia: https://en.wikipedia.org/wiki/Ziggurat_algorithm

For a normal distribution, Marsaglia suggests a compact algorithm:
  1.  Let x = −ln(U1)/x1.
  2.  Let y = −ln(U2).
  3.  If 2y > x^2, return x + x1.
  4.  Otherwise, go back to step 1.
-}
zigguratNormalTailMarsaglia : Float -> Random.Generator Float
zigguratNormalTailMarsaglia x1 =
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

tableSize = 256
normalZigguratTables : Array.Array (Float, Float)
normalZigguratTables =
  let
    n = tableSize
    pdfFunc = normalDensity 0 1
    invPdfFunc = normalDensityInverse 0 1
    cdfFunc = standardNormalCumulative
    x1 = zigguratX1 n pdfFunc invPdfFunc cdfFunc
    y1 = pdfFunc x1
    tailArea = 1 - standardNormalCumulative x1
    layerArea = x1*y1 + tailArea
  in
    Array.fromList <| zigguratTables n y1 layerArea pdfFunc invPdfFunc


makeTwoSided oneSidedGen =
  Random.map2 (\x b -> if b then x else -x) oneSidedGen Random.bool

{-| Generate samples from a standard normal distribution.
-}
normal : Random.Generator Float
-- normal = ziggurat <| normalDensity 0 1
normal =
  let
    tables = normalZigguratTables
    pFunc = normalDensity 0 1
    (x1, y1) = case Array.get 1 tables of
    -- (x1, y1) = case List.head tables of
      Just pair -> pair
      Nothing -> Debug.crash "The ziggurat tables for the normal distribution were empty"
    -- tailFunc = zigguratNormalTail x1
    tailFunc = zigguratNormalTailMarsaglia x1
    oneSidedNormal = ziggurat tables pFunc tailFunc
  in
    makeTwoSided oneSidedNormal


simpleNormal : Random.Generator Float
simpleNormal =
  let
    pFunc = normalDensity 0 1
    invPFunc = normalDensityInverse 0 1
    f0 = pFunc 0
    u1u2gen = Random.pair probability probability
    -- fallback (u1, u2) =
    --   let
    --     xMax = invPFunc (u1 * f0)
    --     x = u2 * xMax
    --   in
    --     identityGenerator x
    fallback (u1, u2) =
      let
        x = probit u1
        -- x = u2 * xMax
      in
        identityGenerator x
    oneSidedNormal = u1u2gen `Random.andThen` fallback
  in
    makeTwoSided oneSidedNormal
