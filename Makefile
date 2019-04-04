compile:
	f2py3 -c boat.f95 -m boat
	f2py3 -c trade.f95 -m trade
	f2py3 -c flight.f95 -m flight

