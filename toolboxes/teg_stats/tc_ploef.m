function tc_ploef(varargin)

% Hallo met Corinde!
%
% Nesting is columns in rows.
% Rows are plotted on the main X axis.
% Columns are plotted within the clusters.
%
% That is: each row will be plotted as cluster.

m = varargin{1};
se = varargin{2};
if length(varargin) > 2,
    legendcell = varargin{3};
else
    legendcell = {}; for n = 1:size(m, 1), legendcell{n} = num2str(n); end;
end;

if size(m, 1) == 1,
    m = m';
    se = se';
end;

xminall = Inf;
xmaxall = -Inf;
yminall = Inf;
ymaxall = -Inf;

b0 = bar(m);
legend(legendcell);
hold on;
for ib = 1:length(b0),
    c0 = get(get(b0(ib), 'Children'));
    tops0 = max(c0.YData);
    bottoms0 = min(c0.YData);
    min0x = min(c0.XData);
    if min(min0x) < xminall,
        xminall = min(min0x);
    end;
    max0x = max(c0.XData);
    if max(max0x) > xmaxall,
        xmaxall = max(max0x);
    end;
    x0v = mean([min0x; max0x]);
    for isb = 1:length(tops0),
        barWidth = max0x(isb) - min0x(isb);
        dw = barWidth * 0.25;
        x0 = x0v(isb);
        if abs(tops0(isb)) > abs(bottoms0(isb)),
            y0 = tops0(isb);
        else,
            y0 = bottoms0(isb);
        end;
        dh = se(isb, ib);
        if abs(tops0(isb)) > abs(bottoms0(isb)),
            plot([x0 x0], [y0 y0 + dh], 'k-');
            plot([x0 - dw x0 + dw], [y0 + dh y0 + dh], 'k-');
        else,
            plot([x0 x0], [y0 y0 - dh], 'k-');
            plot([x0 - dw x0 + dw], [y0 - dh y0 - dh], 'k-');
        end;
        if y0 - dh < yminall,
            yminall = y0 - dh;
        end;
        if y0 + dh > ymaxall,
            ymaxall = y0 + dh;
        end;
    end;
end;

set(gca, 'XTick', 1:size(m, 1));

xbuff = (xmaxall - xminall) * 0.05;
xlim([xminall - xbuff, xmaxall + xbuff]);
ybuff = (ymaxall - yminall) * 0.05;
ylim([yminall - ybuff ymaxall + ybuff]);

fprintf(['print -dtiff -r600 plotfilename\n']);
