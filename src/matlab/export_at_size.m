function export_at_size(fig, name, size, dpi)
	arguments
		fig
		name
		size=[5, 5]
		dpi=300
	end
	set(fig, "Units", "inches");
	set(fig, "Position", [6, 1, size(1), size(2)]);
	exportgraphics(fig, name, resolution=dpi);
end
