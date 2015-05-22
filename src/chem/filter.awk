/.*exposeCando.*/ {
    if ( candoLine != "" ) print candoLine;
    candoLine = $0; }
/.*exposePython.*/ {
    candoLine = "";
}