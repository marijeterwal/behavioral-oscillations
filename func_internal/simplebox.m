function qs = simplebox(labels, dat, colors)
%{ 
----- simplebox -----

DESCRIPTION:
Makes narrow form box plots of groups of data

INPUTS:
- labels: grouping info of data (integers or floats)
- data: data array
- colors: rgb colors to use for the box plots (optional, default: black)

OUTPUTS:
- qs: mean per label

Marije ter Wal - 2020
m.j.terwal@bham.ac.uk

%}

assert(all(size(labels) == size(dat)),'Data and labels do not match')

if isempty(colors)
    colors = [0,0,0];
end

if size(colors,1) == 1
    colors = repmat(colors,[length(unique(labels)),1]);
end

hold on
id = 1;
qs = zeros(length(unique(labels)),1);
for lb = unique(labels, 'stable')'
    datalb = dat(labels==lb);
    
    if length(datalb)>3
        q = quantile(datalb,[0.05 0.25 0.50 0.75 0.95]);
        qs(id) = q(3);

        % plot box
        plot([lb,lb], [q(1),q(5)],'linewidth',1, 'color', colors(id,:)) 
        plot([lb,lb], [q(2),q(4)],'linewidth',3, 'color', colors(id,:)) 
        % plot mean
        scatter(lb, q(3),50, colors(id,:), 'filled') 
    else
        scatter(lb*ones(length(datalb),1), datalb, 30, colors(id,:), 'filled') 
    end
    id = id+1;
end

hold off
end