%hellingerDist measures the hellinger distance between pairs of distributions
%
%
%	Use:
%		takes a reference distribution and a matrix in which each row is used as a 
%		comparison distribution
%
%
% Syntax:
%   [hdists] = hellingerDist(indists, ref, b)
%
% Inputs:
%   indists - a matrix; each row is a distribution to be compared to ref
%   ref 		- an array reference distribution
%		b		 - number of bins.
%
% Outputs:


function [hdists, dists, test] = hellingerDist(indists, ref, b)

	if(nargin < 3)
		error('Too few arguments!\n');
	end

	%set number of bins
	for i=1:length(indists(:,1))%for each row in the matrix
		dist = indists(i,:);
		%find the maximum value
		maxval = 0;
		if (max(ref) > max(dist))
			maxval = max(ref);
		else 
			maxval = max(dist);
		end

		%get minimum value
		minval = 999999;%some big number
		if (min(ref) < min(dist))
			minval = min(ref);
		else
			minval = min(dist);
		end

		%get the bin-width for the separate distributions
		binwidth = maxval / b;

		rbinwidth = max(ref) / b;
		binwidth = max(dist) / b;

		%get bins for test set (linked expression level data) from tdist
		c = 1;
		h = 0;%hellinger distance intermediate variable
		tbins = zeros(1,b);
		binnum = 1;
		for t=1:length(dist)
			pr = (dist(t)/binwidth);
			if (pr == 0)
				pr = 1;
			end
			if (pr > b)%account for experiment3
				pr = b;
			end
			binnum = ceil(pr);
			tbins(binnum) = tbins(binnum) + 1;
		end
			
		%get bins for ref
		c = 1;
		rbins = zeros(1,b);
		for t=1:length(ref)
			pr = int32(ref(t)/rbinwidth);
			
			if (pr == 0)
				pr = 1;
			end
			if (pr > b)%account for experiment3
				pr = b;
			end
			binnum = ceil(pr);
			rbins(binnum) = rbins(binnum) + 1;
		end
			

		%sum the difference between the sqrt of proportion of features in the bins
		%from each population
		h = 0;%hellinger distance intermediate variable
		for n=1:b
			h = h + (sqrt(tbins(n)/length(dist)) - sqrt(rbins(n)/length(ref)))^2;
		end;
		h = sqrt(h);%hellinger distance(xp. lvl. dist at maker i, avg xp. lvl. dist)
		hdists(i) = h;
	end
end
	

