#include "contact/contactKktSystem.h"

#include <cmath>
#include <stdexcept>

namespace
{
void validate(const ContactKktLinearization &linearization)
{
	const Eigen::Index primalDofs = linearization.residual.size();
	const Eigen::Index contacts = linearization.gaps.size();
	if (primalDofs <= 0 ||
		linearization.hessian.rows() != primalDofs ||
		linearization.hessian.cols() != primalDofs ||
		linearization.constraintJacobian.rows() != contacts ||
		linearization.constraintJacobian.cols() != primalDofs ||
		linearization.multipliers.size() != contacts)
	{
		throw std::invalid_argument("Contact KKT linearization dimensions are inconsistent");
	}
	if (!linearization.hessian.allFinite() ||
		!linearization.residual.allFinite() ||
		!linearization.constraintJacobian.allFinite() ||
		!linearization.gaps.allFinite() ||
		!linearization.multipliers.allFinite())
	{
		throw std::invalid_argument("Contact KKT linearization contains nonfinite values");
	}
}
}

ContactKktSystem buildContactKktSystem(
	const ContactKktLinearization &linearization)
{
	validate(linearization);
	const int primalDofs = static_cast<int>(linearization.residual.size());
	const int contacts = static_cast<int>(linearization.gaps.size());

	ContactKktSystem system;
	system.primalDofs = primalDofs;
	system.contactCount = contacts;
	system.matrix = MatrixXd::Zero(
		primalDofs + contacts, primalDofs + contacts);
	system.rightHandSide = VectorXd::Zero(primalDofs + contacts);
	system.matrix.topLeftCorner(primalDofs, primalDofs) =
		linearization.hessian;
	if (contacts > 0)
	{
		system.matrix.topRightCorner(primalDofs, contacts) =
			linearization.constraintJacobian.transpose();
		system.matrix.bottomLeftCorner(contacts, primalDofs) =
			linearization.constraintJacobian;
	}
	system.rightHandSide.head(primalDofs) = -(
		linearization.residual -
		linearization.constraintJacobian.transpose() * linearization.multipliers);
	system.rightHandSide.tail(contacts) = -linearization.gaps;
	return system;
}

ContactKktSolution solveContactKktSystem(
	const ContactKktLinearization &linearization)
{
	const ContactKktSystem system = buildContactKktSystem(linearization);
	ContactKktSolution result;
	Eigen::FullPivLU<MatrixXd> factorization(system.matrix);
	if (!factorization.isInvertible())
	{
		return result;
	}

	const VectorXd solution = factorization.solve(system.rightHandSide);
	if (!solution.allFinite())
	{
		return result;
	}
	result.configurationIncrement = solution.head(system.primalDofs);
	result.multiplierIncrement = -solution.tail(system.contactCount);
	result.updatedMultipliers =
		linearization.multipliers + result.multiplierIncrement;
	result.linearResidualNorm =
		(system.matrix * solution - system.rightHandSide).norm();
	result.success = std::isfinite(result.linearResidualNorm) &&
		result.linearResidualNorm <= 1.0e-9 *
			std::max(1.0, system.rightHandSide.norm());
	return result;
}

ContactActiveSetUpdate updateContactActiveSet(
	const std::vector<ContactCandidate> &activeContacts,
	const VectorXd &multipliers,
	const std::vector<ContactCandidate> &detectedCandidates,
	double gapTolerance,
	double multiplierTolerance)
{
	if (multipliers.size() != static_cast<int>(activeContacts.size()))
	{
		throw std::invalid_argument(
			"Active-contact and multiplier dimensions are inconsistent");
	}
	if (!std::isfinite(gapTolerance) || gapTolerance < 0.0 ||
		!std::isfinite(multiplierTolerance) || multiplierTolerance < 0.0)
	{
		throw std::invalid_argument("Active-set tolerances must be finite and nonnegative");
	}

	ContactActiveSetUpdate result;
	std::vector<double> retainedMultipliers;
	for (int contact = 0; contact < static_cast<int>(activeContacts.size()); ++contact)
	{
		if (multipliers[contact] < -multiplierTolerance)
		{
			++result.released;
			continue;
		}
		result.contacts.push_back(activeContacts[contact]);
		retainedMultipliers.push_back(multipliers[contact]);
	}

	auto contains = [&result](const ContactCandidate &candidate) {
		for (const ContactCandidate &active : result.contacts)
		{
			if (active.rodVertex == candidate.rodVertex &&
				active.boundaryId == candidate.boundaryId)
			{
				return true;
			}
		}
		return false;
	};
	for (const ContactCandidate &candidate : detectedCandidates)
	{
		if (candidate.gap < -gapTolerance && !contains(candidate))
		{
			result.contacts.push_back(candidate);
			retainedMultipliers.push_back(0.0);
			++result.added;
		}
	}
	result.multipliers = VectorXd(static_cast<int>(retainedMultipliers.size()));
	for (int contact = 0; contact < result.multipliers.size(); ++contact)
	{
		result.multipliers[contact] = retainedMultipliers[contact];
	}
	return result;
}
