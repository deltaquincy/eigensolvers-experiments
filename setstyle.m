function setstyle(ax, style)
% Helper function: Set style of a figure.
switch style
    case 'latex'
        ax.FontSize = 15;
        ax.TickLabelInterpreter = 'latex';
        ax.XLabel.Interpreter = 'latex';
        ax.YLabel.Interpreter = 'latex';
        if ~isempty(ax.Title)
            ax.Title.Interpreter = 'latex';
        end
        if ~isempty(ax.Legend)
            ax.Legend.Interpreter = 'latex';
        end
end

