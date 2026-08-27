#include "contact/contactDetector.h"
#include "geometry/planarSlabDomain.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace
{
bool near(double left, double right, double tolerance = 1.0e-12)
{
	return std::abs(left - right) <= tolerance;
}

bool rejectedInvalidParameters(const PlanarSlabDomain &domain)
{
	try
	{
		ContactDetector invalid(domain, -0.01, 0.02, 0.03);
		(void)invalid;
		return false;
	}
	catch (const std::invalid_argument &)
	{
	}
	try
	{
		ContactDetector invalid(domain, 0.01, -0.02, 0.03);
		(void)invalid;
		return false;
	}
	catch (const std::invalid_argument &)
	{
	}
	try
	{
		ContactDetector invalid(domain, 0.01, 0.02, -0.03);
		(void)invalid;
		return false;
	}
	catch (const std::invalid_argument &)
	{
	}
	return true;
}
}

int main()
{
	const PlanarSlabDomain domain(
		Vector3d::Zero(), Vector3d::UnitZ(), 0.1, 0.1);
	const ContactDetector detector(domain, 0.01, 0.02, 0.03);

	const ContactDetectionResult centered = detector.detect(
		{Vector3d(0.0, 0.0, 0.0), Vector3d(0.2, 0.0, 0.0)});
	if (!centered.bodyCandidates.empty() || !near(centered.minimumBodyGap, 0.09) ||
		!centered.tip.safe)
	{
		std::cerr << "Centered-rod contact detection failed.\n";
		return 1;
	}

	const std::vector<Vector3d> vertices = {
		Vector3d(0.0, 0.0, 0.0),
		Vector3d(0.1, 0.0, 0.075),
		Vector3d(0.2, 0.0, -0.085),
		Vector3d(0.3, 0.0, 0.095),
		Vector3d(0.4, 0.0, 0.08)};
	const ContactDetectionResult result = detector.detect(vertices);
	if (result.bodyCandidates.size() != 3 || !near(result.minimumBodyGap, -0.005))
	{
		std::cerr << "Unexpected body-contact candidate count or minimum gap.\n";
		return 1;
	}

	const ContactCandidate &penetrating = result.bodyCandidates[0];
	const ContactCandidate &lower = result.bodyCandidates[1];
	const ContactCandidate &upper = result.bodyCandidates[2];
	if (penetrating.rodVertex != 3 || penetrating.boundaryId != 1 ||
		!near(penetrating.centerlineClearance, 0.005) ||
		!near(penetrating.gap, -0.005) ||
		(penetrating.normal + Vector3d::UnitZ()).norm() > 1.0e-12 ||
		lower.rodVertex != 2 || lower.boundaryId != 0 || !near(lower.gap, 0.005) ||
		(lower.normal - Vector3d::UnitZ()).norm() > 1.0e-12 ||
		upper.rodVertex != 1 || upper.boundaryId != 1 || !near(upper.gap, 0.015) ||
		(upper.normal + Vector3d::UnitZ()).norm() > 1.0e-12)
	{
		std::cerr << "Candidate identity, ordering, gap, or normal failed.\n";
		return 1;
	}

	if (result.tip.safe || !near(result.tip.clearance, 0.02) ||
		!near(result.tip.requiredClearance, 0.03) || result.tip.boundaryId != 1 ||
		(result.tip.normal + Vector3d::UnitZ()).norm() > 1.0e-12)
	{
		std::cerr << "Unsafe tip-clearance detection failed.\n";
		return 1;
	}

	// The last vertex is always the forbidden-contact tip, even when it is
	// deeply inside the body's activation region.
	const ContactDetectionResult tipExcluded = detector.detect(
		{Vector3d(0.0, 0.0, 0.0), Vector3d(0.0, 0.0, 0.099)});
	if (!tipExcluded.bodyCandidates.empty() || tipExcluded.tip.safe)
	{
		std::cerr << "Distal tip was treated as an admissible body contact.\n";
		return 1;
	}

	const ContactDetectionResult singleTip =
		detector.detect({Vector3d(0.0, 0.0, 0.05)});
	if (!singleTip.bodyCandidates.empty() ||
		!std::isinf(singleTip.minimumBodyGap) || !singleTip.tip.safe)
	{
		std::cerr << "Single-tip detection failed.\n";
		return 1;
	}

	bool rejectedEmpty = false;
	try
	{
		detector.detect({});
	}
	catch (const std::invalid_argument &)
	{
		rejectedEmpty = true;
	}
	if (!rejectedEmpty || !rejectedInvalidParameters(domain))
	{
		std::cerr << "Contact detector input validation failed.\n";
		return 1;
	}

	bool rejectedNonfinite = false;
	try
	{
		detector.detect({Vector3d(
			std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0)});
	}
	catch (const std::invalid_argument &)
	{
		rejectedNonfinite = true;
	}
	if (!rejectedNonfinite)
	{
		std::cerr << "Contact detector accepted a nonfinite tip.\n";
		return 1;
	}

	std::cout << "Body-contact candidate and tip-clearance checks passed.\n";
	return 0;
}
