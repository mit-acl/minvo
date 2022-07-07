function plotMatrix(matrix_value, label_colorbar, all_deg, all_seg,font_size)

imagesc(matrix_value); %title (title_string);

%%%% Label Stuff
xlabel('$n$'); ylabel('$s$');
c=colorbar;  c.Label.Interpreter = 'latex'; c.Label.String = label_colorbar; c.TickLabelInterpreter= 'latex'; c.Label.FontSize=13;

caxis([0,2])

tmp=c.TickLabels;
c.TickLabels{end}=['$>',c.TickLabels{end},'$'];

xticklabels(sprintfc('%d',all_deg)); yticklabels(sprintfc('%d',all_seg)); %yticklabels(sprintfc('%d',flip(all_seg)));
set(gca, 'YTick', 1:numel(all_seg)); set(gca, 'XTick', 1:numel(all_deg));
%%%%

% font_size=4.6;

%%%Plot rectangles
hold on;
for i=1:size(matrix_value,1)
    for j=1:size(matrix_value,2)
        
        value_string=num2str(formatNumber(matrix_value(i,j),2));
        if(matrix_value(i,j)>1) %MINVO performs better
            rectangle('Position',[j-0.5,i-0.5,1.0,1.0],'EdgeColor','r','LineWidth',1.3)
%             plot([j-0.5,j+0.5],[i-0.5,i+0.5])
%             plot(j,i,'.r')
            text(j,i,['\textbf{',value_string,'}'],'HorizontalAlignment','center','FontSize',font_size,'Color','r')
        else
            text(j,i,value_string,'HorizontalAlignment','center','FontSize',font_size,'Color','k')
        end
    end
end
end