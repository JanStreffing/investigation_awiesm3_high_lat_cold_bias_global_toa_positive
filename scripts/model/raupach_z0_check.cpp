// Verification of the Raupach (1994) roughness implementation added to LPJ-GUESS.
//
// WHY A SEPARATE CHECK.  The implementation is going into a coupled model where a wrong
// constant would show up only as a temperature bias months later, and the one piece of
// independent evidence we have is the curve in the LPJ-GUESS group's roughness note
// (report/Roughness_length.pdf).  That figure is for a canopy of 20 m: z0 rises to a
// peak near LAI+SAI = 0.6, crosses the IFS curve near 2.5, and falls to ~1.1 at 4.4.
// If this implementation does not reproduce those three points, the constants are wrong.
//
// It also checks the two properties the coupling depends on and that no published curve
// shows: that the drag-space aggregation is not the same as averaging z0, and that it
// is exact when all patches are identical.
//
// Build and run:
//   g++ -O2 -I<lpjg>/framework raupach_z0_check.cpp <lpjg>/framework/raupach_z0.cpp
//       -o /tmp/raupach_check && /tmp/raupach_check

#include "raupach_z0.h"
#include <cstdio>
#include <cmath>

static int failures = 0;

static void expect(const char* what, double got, double want, double tol) {
	const bool ok = std::fabs(got - want) <= tol;
	std::printf("  %-58s %8.4f  (expect %7.4f +- %.3f)  %s\n",
	            what, got, want, tol, ok ? "ok" : "FAIL");
	if (!ok) ++failures;
}

int main() {

	std::printf("\nRaupach (1994) z0 for a 20 m canopy -- against the Lund note's figure\n\n");

	// The three points readable off the published curve.
	expect("LAI+SAI = 0.5  (sparse/leafless)", raupach_z0(20.0, 0.5 * 0.5), 2.47, 0.10);
	expect("LAI+SAI = 0.6  (the peak)",        raupach_z0(20.0, 0.5 * 0.6), 2.57, 0.15);
	expect("LAI+SAI = 4.4  (dense/full leaf)", raupach_z0(20.0, 0.5 * 4.4), 1.07, 0.10);

	std::printf("\nThe sign that matters: z0 must FALL as the canopy fills in\n\n");
	double prev = raupach_z0(20.0, 0.5 * 0.8);
	bool monotone = true;
	for (double pai = 1.0; pai <= 5.0; pai += 0.2) {
		const double z0 = raupach_z0(20.0, 0.5 * pai);
		if (z0 > prev + 1e-12) monotone = false;
		prev = z0;
	}
	std::printf("  %-58s %s\n", "monotone decreasing for LAI+SAI in [0.8, 5.0]",
	            monotone ? "ok" : "FAIL");
	if (!monotone) ++failures;

	// The IFS behaviour we are replacing, for contrast: cvh = 1-exp(-0.5 LAI) scaling a
	// fixed 2.0 m table value.  This is what the model does today.
	std::printf("\nWhat the model does today, for comparison [m]\n\n");
	std::printf("  %8s %12s %12s\n", "LAI+SAI", "Raupach", "IFS(cvh*2.0)");
	for (double pai = 0.5; pai <= 4.5; pai += 0.5) {
		const double ifs = (1.0 - std::exp(-0.5 * pai)) * 2.0;
		std::printf("  %8.1f %12.3f %12.3f\n", pai, raupach_z0(20.0, 0.5 * pai), ifs);
	}

	std::printf("\nA Siberian larch stand, leaf-on vs leaf-off\n");
	std::printf("  0.1 stems/m2, 0.20 m diameter, 15 m tall -> stem lambda = %.3f\n",
	            0.1 * 0.20 * 15.0);
	{
		RaupachAccumulator on, off;
		on.begin_patch();
		on.add_woody(/*lai*/ 2.5, /*dens*/ 0.1, /*diam*/ 0.20, /*h*/ 15.0);
		on.end_patch(1.0);
		off.begin_patch();
		off.add_woody(/*lai*/ 0.0, /*dens*/ 0.1, /*diam*/ 0.20, /*h*/ 15.0);
		off.end_patch(1.0);
		std::printf("  leaf-on  z0 = %.3f m   leaf-off z0 = %.3f m   -> winter is %s\n",
		            on.gridcell_z0(), off.gridcell_z0(),
		            off.gridcell_z0() > on.gridcell_z0() ? "ROUGHER (correct)"
		                                                 : "SMOOTHER (WRONG)");
		if (off.gridcell_z0() <= on.gridcell_z0()) ++failures;
	}

	std::printf("\nDrag-space aggregation\n\n");
	{
		// Identical patches must return exactly the single-patch value.
		RaupachAccumulator uni;
		for (int i = 0; i < 4; ++i) {
			uni.begin_patch();
			uni.add_woody(2.0, 0.1, 0.2, 15.0);
			uni.end_patch(0.25);
		}
		RaupachAccumulator one;
		one.begin_patch();
		one.add_woody(2.0, 0.1, 0.2, 15.0);
		one.end_patch(1.0);
		expect("4 identical patches == 1 patch", uni.gridcell_z0(), one.gridcell_z0(), 1e-9);

		// A rough patch and a smooth one: the drag-space answer must sit BELOW the
		// arithmetic mean of the two roughness lengths, because drag is logarithmic in
		// z0 and the smooth patch pulls harder than a linear average admits.
		RaupachAccumulator mix;
		mix.begin_patch(); mix.add_woody(0.2, 0.10, 0.25, 20.0); mix.end_patch(0.5);
		mix.begin_patch(); mix.add_woody(3.0, 0.01, 0.03,  2.0); mix.end_patch(0.5);

		RaupachAccumulator a, b;
		a.begin_patch(); a.add_woody(0.2, 0.10, 0.25, 20.0); a.end_patch(1.0);
		b.begin_patch(); b.add_woody(3.0, 0.01, 0.03,  2.0); b.end_patch(1.0);
		const double linear = 0.5 * (a.gridcell_z0() + b.gridcell_z0());

		std::printf("  rough patch z0 = %.4f, smooth patch z0 = %.4f\n",
		            a.gridcell_z0(), b.gridcell_z0());
		std::printf("  %-58s %8.4f\n", "arithmetic mean of z0 (WRONG way)", linear);
		std::printf("  %-58s %8.4f\n", "drag-space aggregate (what we send)",
		            mix.gridcell_z0());
		const bool below = mix.gridcell_z0() < linear;
		std::printf("  %-58s %s\n", "drag-space result below the arithmetic mean",
		            below ? "ok" : "FAIL");
		if (!below) ++failures;
	}

	std::printf("\nFallback behaviour\n\n");
	{
		RaupachAccumulator empty;
		empty.begin_patch();
		empty.end_patch(1.0);
		const bool ok = empty.gridcell_z0() == RAUPACH_Z0_MISSING;
		std::printf("  %-58s %s\n", "no woody area -> RAUPACH_Z0_MISSING (IFS falls back)",
		            ok ? "ok" : "FAIL");
		if (!ok) ++failures;

		const bool ok2 = raupach_z0(0.1, 1.0) == RAUPACH_Z0_MISSING;
		std::printf("  %-58s %s\n", "canopy below MIN_HEIGHT -> missing",
		            ok2 ? "ok" : "FAIL");
		if (!ok2) ++failures;
	}

	std::printf("\n%s: %d failure(s)\n\n", failures ? "FAILED" : "PASSED", failures);
	return failures ? 1 : 0;
}
