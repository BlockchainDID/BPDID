// SPDX-License-Identifier: GPL-3.0
pragma solidity >=0.7.0 <0.9.0;

contract ASC{
    struct Attrset{
        string attrset;
        bool used;
    }
    mapping (bytes32 => Attrset) public Attrsets;
    function setStore(string memory attrsets, bytes32 nonce) public returns(bytes32){
        bytes32 SetId = sha256(abi.encodePacked(msg.sender,nonce));
        Attrsets[SetId].attrset = attrsets;
        Attrsets[SetId].used = true;
        return SetId;
    }
    function setQuery(bytes32 SetId) public view returns(bool){
        return Attrsets[SetId].used;
    }
}