function A = grid_graph(w, h)                      
	px = sparse(1:w-1, 2:w, 1, w, w);          
	py = sparse(1:h-1, 2:h, 1, h, h);          
	A = kron(py,speye(w)) + kron(speye(h),px); 
	A = A + A';                                
end