// SPDX-License-Identifier: GPL-3.0
pragma solidity >=0.7.0 <0.9.0;

contract AMC{
    struct params{
        bytes32 gx;
        bytes32 Yapx;
        bytes32 deltax;
        bool glsb;
        bool Yaplsb;
        bool deltalsb;
        bytes32[] URL;
    }
    address CA;
    address owner;
    params public pubparams;

    constructor(address ca) {
        owner = msg.sender;
        CA = ca;
    }
    function changeca(address newaddress) public{
        require(msg.sender == owner);
        CA = newaddress;
    }
    function paramsInit(bytes32 gx, bytes32 Yapx, bytes32 deltax, bool glsb, bool Yaplsb, bool deltalsb) public{
        require(msg.sender == CA);
        pubparams.gx = gx;
        pubparams.Yapx = Yapx;
        pubparams.deltax = deltax;
        pubparams.glsb = glsb;
        pubparams.Yaplsb = Yaplsb;
        pubparams.deltalsb = deltalsb;
    }
    function accountRevoke(bytes32 deltax, bytes32 wu, bool deltalsb) public{
        require(msg.sender == CA);
        pubparams.deltalsb = deltalsb;
        pubparams.deltax = deltax;
        pubparams.URL.push(wu);
    }
    function getURL() public view returns(bytes32[] memory){
        return pubparams.URL;
    }
    function getURLvalue(uint index) public view returns(bytes32){
        return pubparams.URL[index];
    }
}