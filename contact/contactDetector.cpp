#include "contact/contactDetector.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

ContactDetector::ContactDetector(
	const ConfinedDomain &m_domain,
	double m_rodRadius,
	double m_activationDistance,
	double m_tipSafeDistance)
	: domain(&m_domain),
	  rodRadius(m_rodRadius),
	  activationDistance(m_activationDistance),
	  tipSafeDistance(m_tipSafeDistance)
{
	if (!std::isfinite(rodRadius) || rodRadius < 0.0)
	{
		throw std::invalid_argument("Contact detector rod radius must be nonnegative");
	}
	if (!std::isfinite(activationDistance) || activationDistance < 0.0)
	{
		throw std::invalid_argument(
			"Contact detector activation distance must be nonnegative");
	}
	if (!std::isfinite(tipSafeDistance) || tipSafeDistance < 0.0)
	{
		throw std::invalid_argument(
			"Contact detector tip safety distance must be nonnegative");
	}
}

ContactDetectionResult ContactDetector::detect(
	const std::vector<Vector3d> &rodVertices) const
{
	if (rodVertices.empty())
	{
		throw std::invalid_argument("Contact detection requires at least a tip vertex");
	}

	ContactDetectionResult result;
	const int tipIndex = static_cast<int>(rodVertices.size()) - 1;
	for (int vertex = 0; vertex < tipIndex; ++vertex)
	{
		if (!rodVertices[vertex].allFinite())
		{
			throw std::invalid_argument("Contact detection received a nonfinite vertex");
		}
		const DomainQuery query = domain->query(rodVertices[vertex]);
		if (!std::isfinite(query.clearance) || !query.clearanceNormal.allFinite() ||
			!query.clearanceHessian.allFinite() || !query.closestPoint.allFinite())
		{
			throw std::runtime_error("Confined domain returned a nonfinite body query");
		}
		const double gap = query.clearance - rodRadius;
		result.minimumBodyGap = std::min(result.minimumBodyGap, gap);
		if (gap <= activationDistance)
		{
			ContactCandidate candidate;
			candidate.rodVertex = vertex;
			candidate.boundaryId = query.boundaryId;
			candidate.centerlineClearance = query.clearance;
			candidate.gap = gap;
			candidate.rodCenter = rodVertices[vertex];
			candidate.closestBoundaryPoint = query.closestPoint;
			candidate.normal = query.clearanceNormal;
			candidate.gapHessian = query.clearanceHessian;
			result.bodyCandidates.push_back(candidate);
		}
	}

	std::sort(
		result.bodyCandidates.begin(), result.bodyCandidates.end(),
		[](const ContactCandidate &left, const ContactCandidate &right) {
			if (left.gap != right.gap)
			{
				return left.gap < right.gap;
			}
			if (left.rodVertex != right.rodVertex)
			{
				return left.rodVertex < right.rodVertex;
			}
			return left.boundaryId < right.boundaryId;
		});

	if (!rodVertices.back().allFinite())
	{
		throw std::invalid_argument("Contact detection received a nonfinite tip");
	}
	const DomainQuery tipQuery = domain->query(rodVertices.back());
	if (!std::isfinite(tipQuery.clearance) || !tipQuery.clearanceNormal.allFinite() ||
		!tipQuery.clearanceHessian.allFinite() || !tipQuery.closestPoint.allFinite())
	{
		throw std::runtime_error("Confined domain returned a nonfinite tip query");
	}
	result.tip.clearance = tipQuery.clearance;
	result.tip.requiredClearance = tipSafeDistance;
	result.tip.closestBoundaryPoint = tipQuery.closestPoint;
	result.tip.normal = tipQuery.clearanceNormal;
	result.tip.boundaryId = tipQuery.boundaryId;
	result.tip.safe = tipQuery.clearance >= tipSafeDistance;
	return result;
}
