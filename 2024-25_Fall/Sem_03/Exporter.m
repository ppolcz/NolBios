function ret = Exporter(Name)
    mkdir("Results/" + Name)
    ret = @(Fig,ind,p_val) exportgraphics(Fig,sprintf('Results/%s/%03d.jpg',Name,ind));
end
