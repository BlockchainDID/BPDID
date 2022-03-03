// SPDX-License-Identifier: GPL-3.0
pragma solidity >=0.7.0 <0.9.0;

contract ERC{
    struct pkIDP{
        bytes32 Xx;
        bytes32 Yx;
        bool Xlsb;
        bool Ylsb;
    }
    mapping (bytes32 => pkIDP) public IDPlist;
    function entityRegister(bytes32 Xx, bytes32 Yx, bool Xlsb, bool Ylsb) public returns(bytes32){
        bytes32 IDPid = sha256(abi.encodePacked(Xx, Yx, Xlsb, Ylsb));
        IDPlist[IDPid].Xx = Xx;
        IDPlist[IDPid].Yx = Yx;
        IDPlist[IDPid].Xlsb = Xlsb;
        IDPlist[IDPid].Ylsb = Ylsb;
        return IDPid;
    }
}